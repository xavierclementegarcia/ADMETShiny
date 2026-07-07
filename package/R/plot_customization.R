# ---------------------------------------------------------------------------
# plot_customization.R
# Colour-palette helpers for ggplot2 objects.
# ---------------------------------------------------------------------------

#' Colour palette selector UI element
#'
#' A Shiny \code{selectInput} that lets the user choose a colour palette for
#' the plots produced by the application.
#'
#' @param id Character. Input id to use for the select input.
#' @return A \code{shiny.tag} (select input).
#' @export
palette_selector_ui <- function(id) {
  selectInput(
    inputId = id,
    label = "Colour palette",
    choices = list(
      "Default"             = "default",
      "Viridis (continuous)" = "viridis",
      "Magma (continuous)"   = "magma",
      "Inferno (continuous)" = "inferno",
      "Set1 (discrete)"      = "Set1",
      "Set2 (discrete)"      = "Set2",
      "Dark2 (discrete)"     = "Dark2",
      "Spectral (diverging)" = "Spectral"
    ),
    selected = "default"
  )
}

#' Apply a colour palette to a ggplot object
#'
#' Modifies an existing ggplot object by adding a colour/fill scale matching
#' the selected palette, intelligently detecting whether the colour variable
#' is continuous or discrete.
#'
#' @param p A ggplot object.
#' @param palette_name Character. Name of the palette (as returned by
#'   \code{\link{palette_selector_ui}}).
#' @param data The data.frame used to build \code{p}.
#' @param color_col Character. Name of the column used for colouring, or
#'   \code{"None"}.
#' @return A ggplot object (possibly modified).
#' @export
apply_palette <- function(p, palette_name, data, color_col) {

  if (is.null(palette_name) || palette_name == "default" ||
      is.null(color_col) || color_col == "None" || !inherits(p, "ggplot")) {
    return(p)
  }

  if (!color_col %in% names(data)) {
    return(p)
  }

  is_continuous <- is.numeric(data[[color_col]])

  if (palette_name %in% c("viridis", "magma", "inferno")) {
    if (is_continuous) {
      p <- p + ggplot2::scale_color_viridis_c(option = palette_name) +
        ggplot2::scale_fill_viridis_c(option = palette_name)
    } else {
      p <- p + ggplot2::scale_color_viridis_d(option = palette_name) +
        ggplot2::scale_fill_viridis_d(option = palette_name)
    }
  } else if (palette_name %in% c("Set1", "Set2", "Dark2", "Spectral")) {
    if (is_continuous) {
      p <- p + ggplot2::scale_color_distiller(palette = palette_name) +
        ggplot2::scale_fill_distiller(palette = palette_name)
    } else {
      p <- p + ggplot2::scale_color_brewer(palette = palette_name) +
        ggplot2::scale_fill_brewer(palette = palette_name)
    }
  }

  p
}
