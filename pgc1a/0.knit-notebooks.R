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


render_notebook("1.inspect-pgc1a-screen")
render_notebook("2.inspect-gene-compound")
render_notebook("3.combine")
