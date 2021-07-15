output_format <- "github_document"

render_notebook <-
  function(notebook_name, output_suffix = "", ...) {
    output_file <- paste0(notebook_name, output_suffix, ".md")
    
    rmarkdown::render(
      glue::glue("{notebook_name}.Rmd"),
      output_file = output_file,
      output_dir = "knit_notebooks",
      output_format = output_format,
      ...
    )
    
    output_file_rel <- file.path("knit_notebooks", output_file)
    
    read_lines(output_file_rel) %>%
      str_remove_all(file.path(getwd(), "knit_notebooks/")) %>%
      write_lines(output_file_rel)
    
  }


render_notebook("1.inspect-p38a-screen")

render_notebook("2.inspect-single-cell")

system("montage output/distributions/*_conc_1.png -tile 4x -geometry +1+1 output/distributions_conc_1.png")

system("montage output/distributions/*_conc_10.png -tile 4x -geometry +1+1 output/distributions_conc_10.png")

render_notebook("3.figures")
