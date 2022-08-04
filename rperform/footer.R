
# You can modify the PR comment footer here. You can use github markdown e.g.
# emojis like :tada:.
# See `?Rperform::pr_comment`
link <- "https://github.com/analyticalmonk/Rperform#readme"
glue::glue(
  "\n\nFurther explanation regarding interpretation and",
  " methodology can be found in the [documentation]({link})."
)