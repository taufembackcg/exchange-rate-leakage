# ___________________________________________________________________________________
# ______________________ 05_export_results.R (figuras + tabelas) ___________________
# ___________________________________________________________________________________
# NAO depende de 'results_all' estar em memoria: le os .xlsx exportados por
# 04_analysis-v2-testes.R (um arquivo por pais) e monta o results_all combinado.

library(dplyr)
library(tidyr)
library(ggplot2)
library(kableExtra)
library(stringr)
library(readxl)
library(purrr)
library(tikzDevice)

# ─────────────────────────────────────────────────────────────────────────────
# MODO DE EXECUCAO (ponto 22: um .bat por pais + um .bat "all")
#   Rscript 05_export_results.r br    -> so o pais (heatmap, tabelas condensadas
#                                        e tabela completa)
#   Rscript 05_export_results.r all   -> so a tabela entre paises (secao 4,
#                                        precisa dos .xlsx de todos os paises)
#   Rscript 05_export_results.r       -> tudo (comportamento anterior)
# ─────────────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", ca, value = TRUE)
  if (length(hit) > 0) return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
  getwd()
}
setwd(get_script_dir())
# Variante da rodada (M1 de reviews_v7.md; RUN_VARIANT no .bat): resultados em
# outputs/ (main) ou outputs_<variante>/.
OUT <- local({ v <- tolower(trimws(Sys.getenv("RUN_VARIANT", "main")))
               if (!nzchar(v) || v == "main") "outputs" else paste0("outputs_", v) })

# Paises: os que tem outputs/Resultados_Todos_<pais>.xlsx (ponto 36: o mesmo
# script serve ao conjunto BRICS e ao conjunto de Colombo e Pelagatti, cada um
# na sua pasta de execucao). Sem arquivos, o conjunto BRICS.
all_countries <- sub("^Resultados_Todos_(.*)[.]xlsx$", "\\1",
                     list.files(OUT, pattern = "^Resultados_Todos_.*[.]xlsx$"))
if (length(all_countries) == 0) all_countries <- c("br", "ru", "in", "cn", "za", "id")
run_args <- commandArgs(trailingOnly = TRUE)
run_mode <- if (length(run_args) >= 1 && nzchar(run_args[1])) tolower(trimws(run_args[1])) else "full"
if (!run_mode %in% c("full", "all", all_countries)) {
  stop("Argumento invalido: ", run_mode, ". Use um pais (", paste(all_countries, collapse = ", "),
       "), 'all' ou nenhum.", call. = FALSE)
}
do_country <- run_mode != "all"              # secoes 1-3 (por pais)
do_cross   <- run_mode %in% c("full", "all")  # secao 4 (entre paises)
cat("05_export_results | modo:", run_mode, "\n")

dir.create(file.path(OUT, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "figures_tikz"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "tables"),  showWarnings = FALSE, recursive = TRUE)

# ___________________________________________________________________________________
# _____________________ 0) MONTA results_all A PARTIR DOS XLSX ______________________
# ___________________________________________________________________________________
# Pais unico (modo <pais>) ou todos (modos "all" e "full")
country_list <- if (run_mode %in% all_countries) run_mode else all_countries

read_country_results <- function(cty) {
  fname <- paste0(OUT, "/Resultados_Todos_", cty, ".xlsx")
  if (!file.exists(fname)) {
    warning(paste0("Arquivo nao encontrado: ", fname, " -- pais ignorado."))
    return(NULL)
  }
  read_excel(fname, sheet = "Resultados_Todos")
}

results_all <- map_dfr(country_list, read_country_results)

if (nrow(results_all) == 0) {
  stop("Nenhum resultado encontrado. Verifique se os arquivos outputs/Resultados_Todos_<cty>.xlsx existem.")
}

# ___________________________________________________________________________________
# _____________________________ PARAMETROS DE CONTROLE ______________________________
# ___________________________________________________________________________________
horizons_body <- c("h1", "h4", "h8", "h12")
# Variantes de vazamento (M1) rodam so os horizontes de ONLY_HORIZONS
# (h1, h3, h6, h12): a tabela do corpo usa esses horizontes.
only_h <- trimws(strsplit(Sys.getenv("ONLY_HORIZONS", ""), ",")[[1]])
if (length(only_h)) horizons_body <- only_h

# Nomes dos modelos (problema 13): nomes de Colombo & Pelagatti (2020, Secao
# 4.1), na ordem da Tabela 1; o modelo de precos flexiveis com juros
# (Frenkel-Bilson, ausente em Colombo) recebe nome no mesmo estilo.
model_labels <- c(
  monetary_standard      = "Monetary model",
  monetary_sticky        = "Monetary model, sticky prices",
  monetary_sticky_uip    = "Monetary model, sticky prices, UIP deviations",
  monetary_flexible_uip  = "Monetary model, UIP deviations",
  taylor                 = "Taylor rule model"
)

country_names <- c(
  br = "Brazil",
  cn = "China",
  id = "Indonesia",
  `in` = "India",
  ru = "Russia",
  za = "South Africa",
  # conjunto de Colombo e Pelagatti (ponto 36)
  ca = "Canada",
  jp = "Japan",
  gb = "United Kingdom",
  fr = "France",
  dk = "Denmark",   # ponto 39
  se = "Sweden"
)

sig_stars <- function(p) {
  case_when(
    is.na(p)   ~ "",
    p < 0.01   ~ "***",
    p < 0.05   ~ "**",
    p < 0.10   ~ "*",
    TRUE       ~ ""
  )
}

# ___________________________________________________________________________________
# ______________________________ 1) FIGURA: HEATMAP __________________________________
# ___________________________________________________________________________________

plot_dm_heatmap <- function(df, cty, horizons = NULL) {
  
  d <- df %>%
    filter(country == cty) %>%
    { if (!is.null(horizons)) filter(., horizon %in% horizons) else . }
  
  if (nrow(d) == 0) return(invisible(NULL))
  
  # Seleciona o melhor resultado (menor RMSE) por modelo x horizonte,
  # buscando entre TODAS as combinacoes de method e model_type (dif/ecm)
  d <- d %>%
    group_by(model, horizon) %>%
    slice_min(rmse, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  d <- d %>%
    mutate(
      model   = factor(model_labels[model], levels = model_labels),
      horizon = factor(horizon, levels = paste0("h", 1:12)),
      sig_cat = case_when(
        is.na(dm_test_p_value)   ~ "N/A",
        dm_test_p_value < 0.01   ~ "p < 0.01",
        dm_test_p_value < 0.05   ~ "p < 0.05",
        dm_test_p_value < 0.10   ~ "p < 0.10",
        TRUE                     ~ "Not significant"
      ),
      sig_cat = factor(sig_cat, levels = c("p < 0.01", "p < 0.05", "p < 0.10",
                                           "Not significant", "N/A")),
      # Rotulo da celula: metodo vencedor + (dif/ecm) em nova linha
      label_method = paste0(method, "\n(", model_type, ")")
    )
  
  p <- ggplot(d, aes(x = horizon, y = model, fill = sig_cat)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label_method), size = 2.6, color = "black", fontface = "bold",
              lineheight = 0.85) +
    scale_fill_manual(
      values = c("p < 0.01" = "#08519c", "p < 0.05" = "#3182bd",
                 "p < 0.10" = "#9ecae1", "Not significant" = "#f0f0f0",
                 "N/A" = "grey85"),
      name = "DM test significance"
    ) +
    labs(
      # CORRIGIDO: titulo agora usa nome completo do pais + codigo; subtitulo
      # removido da figura (o texto foi movido para a legenda/caption no .tex).
      title = paste0("DM test significance - ", country_names[cty], " (", toupper(cty), ")"),
      x = "Forecast horizon", y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      # CORRIGIDO: legenda movida para baixo, em uma unica linha horizontal,
      # para que o painel do heatmap use toda a largura disponivel da figura.
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

  fname <- paste0(OUT, "/figures/heatmap_dm_", cty, ".png")
  # CORRIGIDO: largura aumentada e altura reduzida, refletindo o espaco
  # horizontal liberado pela legenda (antes a direita, agora abaixo).
  ggsave(fname, p, width = 9, height = 3.6, dpi = 300)

  # ─────────────────────────────────────────────
  # VERSAO TIKZ (teste solicitado): mesma figura, exportada como codigo TikZ
  # nativo em vez de PNG. Correcoes em relacao ao teste inicial:
  #
  # 1) O bug do "p < 0.01" virando "p ¡ 0.01" NAO era do tikzDevice -- era o
  #    documento de exibicao usando a codificacao OT1 (padrao do LaTeX), que
  #    nao tem o glifo "<". Corrigido no documento com \usepackage[T1]{fontenc}
  #    (o artigo principal ja carrega T1; sanitize = TRUE aqui e mantido como
  #    protecao extra para outros caracteres especiais, mas nao e o que
  #    resolve o "<").
  #
  # 2) Gerado no MESMO tamanho absoluto usado para o PNG (9 x 3.6 in), e nao
  #    no tamanho final de exibicao. Tentativa de gerar direto em 5.1in
  #    (~\textwidth do artigo) fez o texto de cada celula colidir: no TikZ o
  #    texto e vetorial, desenhado em unidades absolutas (pt/mm) que NAO
  #    encolhem com o tamanho do canvas -- reduzir o canvas sem reduzir a
  #    fonte apenas aumenta a fonte em proporcao ao espaco disponivel. Gerando
  #    no tamanho ja validado visualmente pela versao PNG, e encolhendo para
  #    a largura final com \resizebox no documento de exibicao (que escala
  #    tudo, inclusive o texto, proporcionalmente), preserva-se a mesma
  #    relacao texto/celula que ja funcionou no PNG.
  tikz_fname <- paste0(OUT, "/figures_tikz/heatmap_dm_", cty, ".tex")
  tikz(
    file = tikz_fname,
    width = 9, height = 3.6,
    standAlone = FALSE,
    sanitize = TRUE
  )
  print(p)
  dev.off()

  invisible(p)
}

countries <- unique(results_all$country)
if (do_country) walk(countries, ~ plot_dm_heatmap(results_all, .x))

combos <- results_all %>% distinct(country, model_type)

# ___________________________________________________________________________________
# ______________________ 2) TABELA CONDENSADA (CORPO DO TEXTO) ______________________
# ___________________________________________________________________________________

make_condensed_table <- function(df, cty, mtype, horizons = horizons_body) {
  
  tab <- df %>%
    filter(country == cty, model_type == mtype, horizon %in% horizons) %>%
    group_by(model, horizon) %>%
    slice_min(rmse, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  if (nrow(tab) == 0) return(invisible(NULL))
  
  tab <- tab %>%
    mutate(
      cell = sprintf("%.3f%s", relative_rmse, sig_stars(dm_test_p_value)),
      model = model_labels[model],
      horizon = factor(horizon, levels = horizons)
    ) %>%
    select(model, horizon, cell) %>%
    pivot_wider(names_from = horizon, values_from = cell) %>%
    # CORRIGIDO: pivot_wider nao respeitava a ordem dos niveis do fator
    # 'horizon' e ordenava as colunas alfabeticamente (h1, h12, h4, h8, ...);
    # forca explicitamente a ordem correta (h1, h4, h8, h12) via select().
    select(model, all_of(horizons)) %>%
    arrange(match(model, model_labels))

  colnames(tab)[1] <- "Model"
  
  # CORRIGIDO: longtable = TRUE, para ser consistente com make_full_table() e
  # make_mcs_table(). kableExtra emite notas de rodape de duas formas
  # incompativeis dependendo de longtable: com longtable = FALSE, usa
  # \begin{threeparttable} (pacote threeparttable); com longtable = TRUE, usa
  # \begin{ThreePartTable} (pacote threeparttablex). Misturar os dois estilos
  # no mesmo documento LaTeX causa erro fatal ("Illegal nested threeparttable
  # environments"). Como as tabelas completas (make_full_table) e do MCS
  # (make_mcs_table) ja usam longtable = TRUE, esta tabela precisa do mesmo
  # estilo para poder ser \input{} no mesmo artigo.
  tex <- kbl(
    tab, format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
    caption = paste0("Relative RMSE and Diebold-Mariano test significance -- ",
                     toupper(cty), " (", mtype, ")"),
    # CORRIGIDO: kableExtra ja prefixa "tab:" internamente ao label; manter o
    # prefixo aqui gerava \label{tab:tab:dm_...} no .tex final.
    label = paste0("dm_", cty, "_", mtype),
    align = c("l", rep("c", ncol(tab) - 1))
  ) %>%
    kable_styling(latex_options = c("repeat_header")) %>%
    footnote(
      general = "Values are relative RMSE (model / random walk). *, **, *** denote significance at 10%, 5%, 1% under the DM test (H0: RW is not outperformed).",
      general_title = "Note:", footnote_as_chunk = TRUE, threeparttable = TRUE
    )
  
  fname <- paste0(OUT, "/tables/table_condensed_", cty, "_", mtype, ".tex")
  writeLines(as.character(tex), fname)
  invisible(tex)
}

if (do_country) pwalk(combos, ~ make_condensed_table(results_all, ..1, ..2))

# ___________________________________________________________________________________
# ________________________ 3) TABELA COMPLETA (APENDICE) ____________________________
# ___________________________________________________________________________________

# Codigos curtos dos modelos na tabela completa (ponto 35 de review_problems.md):
# com os nomes de Colombo e a coluna CW, a coluna Model tinha de ter 2.2cm e os
# nomes quebravam em tres linhas (Material Suplementar com 518 paginas). A
# legenda vai na nota da tabela.
model_codes <- c(
  monetary_standard      = "M",
  monetary_sticky        = "M-SP",
  monetary_sticky_uip    = "M-SP-UIP",
  monetary_flexible_uip  = "M-UIP",
  taylor                 = "TR"
)
model_codes_legend <- paste(
  "Model codes: M = monetary model; M-SP = monetary model, sticky prices;",
  "M-SP-UIP = monetary model, sticky prices, UIP deviations; M-UIP = monetary model, UIP deviations;",
  "TR = Taylor rule model."
)

make_full_table <- function(df, cty) {
  if (!"cw_test_p_value" %in% names(df)) df$cw_test_p_value <- NA_real_
  
  tab <- df %>%
    filter(country == cty) %>%
    # CORRIGIDO: horizon era ordenado como texto (h1, h10, h11, h12, h2, h3, ...);
    # convertido para fator com niveis numericos antes do arrange/select.
    mutate(horizon = factor(horizon, levels = paste0("h", 1:12))) %>%
    arrange(model_type, model, horizon, method) %>%
    mutate(
      model   = unname(model_codes[model]),
      rmse    = sprintf("%.4f", rmse),
      relative_rmse = sprintf("%.3f", relative_rmse),
      dm_test_p_value = sprintf("%.3f%s", dm_test_p_value, sig_stars(dm_test_p_value)),
      # ponto 30: Clark-West (modelos aninhados), ao lado do DM
      cw_test_p_value = sprintf("%.3f%s", cw_test_p_value, sig_stars(cw_test_p_value))
    ) %>%
    select(
      `Model type` = model_type, Model = model, Horizon = horizon,
      Method = method, RMSE = rmse, `Rel. RMSE` = relative_rmse,
      `DM p-value` = dm_test_p_value, `CW p-value` = cw_test_p_value
    )
  
  # CORRIGIDO: o linesep padrao do kableExtra insere \addlinespace a cada 5
  # linhas, sem relacao com o agrupamento real dos dados (cada combinacao de
  # Model type x Model x Horizon tem 13 linhas, uma por metodo). Aqui o
  # espaco extra e inserido apenas quando esse agrupamento muda de fato.
  group_key <- paste(tab$`Model type`, tab$Model, tab$Horizon)
  n_tab <- nrow(tab)
  linesep_vec <- ifelse(group_key[-n_tab] != group_key[-1], "\\addlinespace", "")

  tex <- kbl(
    tab, format = "latex", booktabs = TRUE, longtable = TRUE, escape = FALSE,
    caption = paste0("Full results (RMSE, relative RMSE, DM and CW tests) -- ", toupper(cty)),
    # CORRIGIDO: mesmo motivo do label acima (evita \label{tab:tab:full_...}).
    label = paste0("full_", cty),
    # Ponto 35: codigos curtos (model_codes) cabem numa linha; espacamento
    # entre colunas de 3pt (abaixo).
    align = c("l", "l", "c", "l", "c", "c", "c", "c"),
    linesep = linesep_vec
  ) %>%
    kable_styling(latex_options = c("repeat_header"), font_size = 8) %>%
    footnote(
      general = paste("*, **, *** denote significance at 10%, 5%, 1% under the Diebold-Mariano (DM) and",
                      "Clark-West (CW) tests (H0: RW is not outperformed; one-sided).", model_codes_legend),
      general_title = "Note:", footnote_as_chunk = TRUE, threeparttable = TRUE
    )

  tex <- sub("\\selectfont", "\\selectfont\\setlength{\\tabcolsep}{3pt}", as.character(tex), fixed = TRUE)
  fname <- paste0(OUT, "/tables/table_full_", cty, ".tex")
  writeLines(as.character(tex), fname)
  invisible(tex)
}

if (do_country) walk(unique(results_all$country), ~ make_full_table(results_all, .x))

# ___________________________________________________________________________________
# _________ 4) MAGNITUDE: RMSE RELATIVA POR PAIS E HORIZONTE (problema 8) ___________
# ___________________________________________________________________________________
# Problema 8 de review_problems.md: as taxas de vitoria (% de especificacoes com
# RMSE relativa < 1) sao descritivas e nao mostram o TAMANHO do ganho. Esta
# secao acrescenta (sem remover nenhuma tabela existente) a distribuicao da RMSE
# relativa entre todas as especificacoes (modelo economico x especificacao x
# metodo), por pais e horizonte, na janela comum de avaliacao (problema 0).
#
# Saidas:
#   outputs/tables/relrmse_country_horizon.csv  -- todas as estatisticas, mais
#       uma linha "all" por pais (todos os horizontes juntos)
#   outputs/tables/table_relrmse_country_horizon.tex -- tabela no estilo do
#       artigo: mediana e, abaixo, intervalo interquartil entre colchetes.

relrmse_stats <- function(d) {
  d %>%
    summarise(
      n_specs      = n(),
      n_obs        = if ("n_obs" %in% names(d)) max(n_obs) else NA_real_,
      mean         = mean(relative_rmse),
      sd           = sd(relative_rmse),
      min          = min(relative_rmse),
      p25          = unname(quantile(relative_rmse, 0.25)),
      median       = median(relative_rmse),
      p75          = unname(quantile(relative_rmse, 0.75)),
      max          = max(relative_rmse),
      pct_below_1  = 100 * mean(relative_rmse < 1),
      pct_dm_5     = 100 * mean(dm_test_p_value < 0.05, na.rm = TRUE),
      pct_cw_5     = if ("cw_test_p_value" %in% names(d)) 100 * mean(cw_test_p_value < 0.05, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
}

if (do_cross) {

rr <- results_all %>% filter(is.finite(relative_rmse))

rr_by_h <- rr %>%
  group_by(country, horizon) %>%
  relrmse_stats() %>%
  mutate(horizon = factor(horizon, levels = paste0("h", 1:12)))

rr_all_h <- rr %>%
  group_by(country) %>%
  relrmse_stats() %>%
  mutate(horizon = "all", n_obs = NA_real_)

rr_csv <- bind_rows(rr_by_h %>% mutate(horizon = as.character(horizon)), rr_all_h) %>%
  mutate(country_name = unname(country_names[country]),
         horizon = factor(horizon, levels = c(paste0("h", 1:12), "all"))) %>%
  arrange(country_name, horizon) %>%
  select(country, country_name, horizon, everything())

write.csv(rr_csv, file.path(OUT, "tables/relrmse_country_horizon.csv"), row.names = FALSE)

make_relrmse_table <- function(st) {
  ctys <- st %>% distinct(country) %>%
    mutate(name = unname(country_names[country])) %>% arrange(name)

  cell_med <- function(cty, h) {
    v <- st$median[st$country == cty & st$horizon == h]
    if (length(v) == 0) "--" else sprintf("%.3f", v)
  }
  cell_iqr <- function(cty, h) {
    i <- which(st$country == cty & st$horizon == h)
    if (length(i) == 0) "" else sprintf("[%.2f, %.2f]", st$p25[i], st$p75[i])
  }

  block <- function(hs) {
    head <- paste0("Country & ", paste(hs, collapse = " & "), " \\\\")
    rows <- unlist(lapply(seq_len(nrow(ctys)), function(k) {
      cty <- ctys$country[k]
      c(paste0(ctys$name[k], " & ", paste(sapply(hs, cell_med, cty = cty), collapse = " & "), " \\\\"),
        paste0(" & ", paste(sapply(hs, cell_iqr, cty = cty), collapse = " & "), " \\\\"),
        if (k < nrow(ctys)) "\\addlinespace" else NULL)
    }))
    c("\\toprule", head, "\\midrule", rows, "\\bottomrule")
  }

  c("\\begin{table*}[htbp]",
    "\\centering",
    "\\caption{Relative RMSE by country and horizon}",
    "\\label{tab:relrmse_country_horizon}",
    "{\\fontsize{8}{10}\\selectfont\\setlength{\\tabcolsep}{3pt}",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lcccccc}",
    block(paste0("h", 1:6)),
    "\\end{tabular*}",
    "",
    "\\vspace{6pt}",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lcccccc}",
    block(paste0("h", 7:12)),
    "\\end{tabular*}",
    "}",
    "\\par\\noindent\\begin{minipage}{\\textwidth}",
    "{\\fontsize{8}{10}\\selectfont Note: median relative RMSE (model RMSE divided by the random walk RMSE) across all specifications (economic model, specification and statistical learning method), evaluated over the common forecast window of each country and horizon. The interquartile range is in brackets. Values below one indicate a smaller RMSE than the random walk. The Taylor model is not estimated for Indonesia.}",
    "\\end{minipage}",
    "\\end{table*}")
}

writeLines(make_relrmse_table(rr_by_h), file.path(OUT, "tables/table_relrmse_country_horizon.tex"))

}  # fim do if (do_cross)

cat("Figuras salvas em outputs/figures/ | Tabelas salvas em outputs/tables/\n")