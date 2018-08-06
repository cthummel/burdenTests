FILE(REMOVE_RECURSE
  "CMakeFiles/gsl_project"
  "CMakeFiles/gsl_project-complete"
  "gsl/src/gsl_project-stamp/gsl_project-install"
  "gsl/src/gsl_project-stamp/gsl_project-mkdir"
  "gsl/src/gsl_project-stamp/gsl_project-download"
  "gsl/src/gsl_project-stamp/gsl_project-update"
  "gsl/src/gsl_project-stamp/gsl_project-patch"
  "gsl/src/gsl_project-stamp/gsl_project-configure"
  "gsl/src/gsl_project-stamp/gsl_project-build"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/gsl_project.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
