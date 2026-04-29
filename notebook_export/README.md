# notebook_export/

This directory contains a pre-rendered HTML export of `PROJECT_STORY.Rmd`.

| File | Description |
|---|---|
| `PROJECT_STORY.html` | Self-contained HTML rendering of the project narrative. Open in any browser — no R installation required. |

## How to regenerate

```r
# Requires R ≥ 4.2 with rmarkdown and all packages in docs/environment.yml
rmarkdown::render(
  "PROJECT_STORY.Rmd",
  output_format = "html_document",
  output_dir    = "notebook_export"
)
```

The HTML file is not tracked in this repository because it is a build artifact and can exceed GitHub's file-size recommendations. Generate it locally using the command above.
