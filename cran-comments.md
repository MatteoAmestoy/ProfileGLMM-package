## Resubmission summary
This is a resubmission. Significant updates include:

* **S3 Integration**: The main functions now returns formal S3 objects of 
  class 'pglmm_data', 'pglmm_mcmc' and 'pglmm_fit'.
* **New Methods**: Added print S3 methods for all the objects.
Added a 'predict', 'plot', and 'summary' S3 methods for 'pglmm_fit' objects to improve 
the user workflow.
* **New Documentation**:
Added a new Introduction Vignette (Intro_to_ProfileGLMM) to provide a clear entry point for new users, covering the basic workflow and S3 method usage.
Updated all .Rd files to reflect the new S3 structure.
