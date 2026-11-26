# spline_glmnet(): modelo caret de SP/SPI (Colombo e Pelagatti, 2020, Secao 4.2.1).
# Usado por 03_models.R e 03c_splines.R (arquivo unico para os dois ficarem iguais).
#
# ─────────────────────────────────────────────
# 2) B-SPLINES COM NOS POR JANELA (Secao 4.2.1 de Colombo & Pelagatti, 2020)
# ─────────────────────────────────────────────
# CORRECAO: problema 19 de review_problems.md. Antes a base B-spline era
# construida uma unica vez sobre o painel inteiro (build_spline_data): os nos
# internos (quantis) e os nos de fronteira (amplitude) usavam a distribuicao
# dos preditores apos cada origem de previsao.
#
# Agora SP e SPI sao um modelo caret proprio (spline_glmnet). Em cada janela:
#   1) para cada preditor, nos internos nos quantis (K = 4, igualmente
#      espacados) e nos de fronteira na amplitude DAQUELA janela; preditor
#      constante na janela (ex.: cambio fixo) nao gera base;
#   2) base cubica aditiva (SP) ou aditiva + produtos entre pares (SPI);
#   3) padronizacao das colunas com media e desvio da janela;
#   4) elastic net (glmnet) com alpha e lambda fixos (tuning na 1a janela).
# A previsao reaplica a base e a padronizacao da janela ao ponto previsto;
# ponto 34 de review_problems.md: valores do ponto previsto fora da amplitude
# da janela sao truncados nos nos de fronteira, e cada coluna da base (inclusive
# os produtos do SPI) e truncada na sua amplitude na janela.
# As funcoes auxiliares ficam DENTRO de spline_glmnet() para viajarem com o
# modelo ate os processos paralelos (o cluster PSOCK nao ve o ambiente global).
# Dummies de quebra do M1 (M6 de reviews_v7.md): colunas so com 0 e 1 entram
# lineares, sem base B-spline (nao ha quantis internos) e fora dos produtos do SPI.
# global_x (variante "scaling", M1 de reviews_v7.md): se informado, os nos, as
# amplitudes e a padronizacao da base vem da amostra inteira (global_x), e nao
# da janela -- reproduz o vazamento; no main fica NULL.
spline_glmnet <- function(interactions = FALSE, K = 4, global_x = NULL) {
# Avalia os argumentos ja aqui: o modelo viaja para os processos paralelos e uma
# promessa nao avaliada (ex.: global_x = sp_global) seria procurada no ambiente
# do processo filho, onde o objeto nao existe.
force(interactions); force(K); force(global_x)

spline_basis_spec <- function(x, K = 4) {
  if (all(x %in% c(0, 1))) return(list(linear = TRUE))       # dummy: entra linear
  lo <- min(x); hi <- max(x)
  if (!is.finite(lo) || hi <= lo) return(NULL)               # constante na janela
  q <- stats::quantile(x, probs = seq(0, 1, length.out = K + 2)[-c(1, K + 2)], names = FALSE)
  list(knots = unique(q[q > lo & q < hi]), boundary = c(lo, hi))
}

spline_design <- function(x, specs, interactions) {
  x <- as.data.frame(x)
  bl <- list()
  lin <- list()
  for (v in names(specs)) {
    s <- specs[[v]]
    if (is.null(s)) next
    if (isTRUE(s$linear)) {
      lin[[v]] <- matrix(as.numeric(x[[v]]), ncol = 1, dimnames = list(NULL, v))
      next
    }
    # Ponto 34: valores fora da amplitude da janela sao truncados nos nos de
    # fronteira (extrapolacao constante). Sem isso a base cubica cresce como
    # polinomio fora da janela e as previsoes de SP/SPI explodiam.
    xv <- pmin(pmax(x[[v]], s$boundary[1]), s$boundary[2])
    b <- suppressWarnings(splines::bs(xv, knots = s$knots, Boundary.knots = s$boundary,
                                      degree = 3, intercept = FALSE))
    b <- unclass(b)[, , drop = FALSE]
    attributes(b) <- list(dim = dim(b))
    colnames(b) <- paste0(v, "_bs", seq_len(ncol(b)))
    bl[[v]] <- b
  }
  X <- do.call(cbind, bl)
  if (interactions && length(bl) >= 2) {
    for (p in combn(names(bl), 2, simplify = FALSE)) {
      B1 <- bl[[p[1]]]; B2 <- bl[[p[2]]]
      for (j in seq_len(ncol(B1))) for (k in seq_len(ncol(B2))) {
        X <- cbind(X, B1[, j] * B2[, k])
        colnames(X)[ncol(X)] <- paste0(p[1], "_bs", j, "X", p[2], "_bs", k)
      }
    }
  }
  if (length(lin) > 0) X <- cbind(X, do.call(cbind, lin))
  X
}

  list(
    label      = if (interactions) "Spline (interactions) + elastic net" else "Spline + elastic net",
    library    = c("glmnet", "splines"),
    type       = "Regression",
    parameters = data.frame(parameter = c("alpha", "lambda"), class = c("numeric", "numeric"),
                            label = c("Mixing Percentage", "Regularization Parameter")),
    grid = function(x, y, len = NULL, search = "grid")
      expand.grid(alpha = seq(0, 1, length.out = 11), lambda = 10^seq(-4, 1, length.out = 50)),
    # Como no glmnet do caret: um ajuste por alpha (caminho inteiro de lambda);
    # os demais lambdas da grade sao "submodelos" previstos do mesmo ajuste.
    loop = function(grid) {
      alphas    <- unique(grid$alpha)
      loop      <- data.frame(alpha = alphas,
                              lambda = sapply(alphas, function(a) max(grid$lambda[grid$alpha == a])))
      submodels <- lapply(alphas, function(a) {
        l <- grid$lambda[grid$alpha == a]
        data.frame(lambda = l[l != max(l)])
      })
      list(loop = loop, submodels = submodels)
    },
    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
      x     <- as.data.frame(x)
      # Nos e padronizacao: da janela (main) ou da amostra inteira (variante scaling)
      xs    <- if (is.null(global_x)) x else as.data.frame(global_x)[, names(x), drop = FALSE]
      specs <- lapply(xs, spline_basis_spec, K = K)
      # dummy constante na janela nao entra (sem variacao para estimar)
      for (v in names(specs)) if (isTRUE(specs[[v]]$linear) && length(unique(x[[v]])) < 2) specs[v] <- list(NULL)
      X     <- spline_design(x, specs, interactions)
      Xs    <- if (is.null(global_x)) X else spline_design(xs, specs, interactions)
      mu    <- colMeans(Xs)
      sdv   <- apply(Xs, 2, stats::sd)
      keep  <- is.finite(sdv) & sdv > 0 & apply(X, 2, stats::sd) > 0
      # Ponto 34: amplitude de cada coluna da base na janela (usada na previsao)
      cmin  <- apply(Xs[, keep, drop = FALSE], 2, min)
      cmax  <- apply(Xs[, keep, drop = FALSE], 2, max)
      X     <- scale(X[, keep, drop = FALSE], center = mu[keep], scale = sdv[keep])
      # Variante scaling: a padronizacao e a da amostra inteira (acima); o glmnet nao
      # repadroniza com a janela. Antes ele repadronizava: coluna da base quase nula
      # na janela (suporte fora dela) ganhava coeficiente enorme e a previsao
      # explodia (relRMSE ate 2e7). No main o X ja vem padronizado pela janela e
      # standardize = TRUE nao muda nada (mesmo resultado de antes).
      fit   <- glmnet::glmnet(X, y, alpha = param$alpha, family = "gaussian",
                              standardize = is.null(global_x))
      list(fit = fit, lambda = param$lambda, specs = specs, keep = keep,
           mu = mu[keep], sdv = sdv[keep], cmin = cmin, cmax = cmax,
           interactions = interactions)
    },
    predict = function(modelFit, newdata, submodels = NULL) {
      X <- spline_design(newdata, modelFit$specs, modelFit$interactions)
      X <- X[, modelFit$keep, drop = FALSE]
      # Ponto 34: cada coluna da base truncada na sua amplitude na janela. Para SP
      # nao muda nada (os preditores ja foram truncados); para SPI evita que um
      # produto quase nulo na janela (desvio minusculo) gere um valor padronizado
      # enorme no ponto previsto.
      if (!is.null(modelFit$cmin)) {
        X <- sweep(sweep(X, 2, modelFit$cmin, pmax), 2, modelFit$cmax, pmin)
      }
      X <- scale(X, center = modelFit$mu, scale = modelFit$sdv)
      out <- as.numeric(stats::predict(modelFit$fit, newx = X, s = modelFit$lambda))
      if (!is.null(submodels) && nrow(submodels) > 0) {
        tmp <- stats::predict(modelFit$fit, newx = X, s = submodels$lambda)
        out <- c(list(out), lapply(seq_len(ncol(tmp)), function(j) as.numeric(tmp[, j])))
      }
      out
    },
    prob = NULL,
    sort = function(x) x[order(-x$lambda, x$alpha), ]
  )
}
