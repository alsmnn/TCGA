# TCGA

Welcome to my analysis project of [`The Cancer Genome Atlas`](https://portal.gdc.cancer.gov/){target="_blank"}!

This little project makes use of [`ProjectTemplate`](http://projecttemplate.net/index.html){target="_blank"}
in order to automate the boring stuff of any project:

* it creates the directory structure and populates it with some information on 
what to do with each directory
* it manages all your needed packages in a central location `config/global.dcf`
* it can load and preprocess your data, enabling you to start with the important
stuff - the analyses
* you can cache intermediate files or big data structures to minimize reload times

Following some contents of the README.md which would be created by
`ProjectTemplate`:

> ProjectTemplate is an R package that helps you organize your statistical
analysis projects. Since you're reading this file, we'll assume that you've
already called `create.project()` to set up this project and all of its
contents.
To load your new project, you'll first need to `setwd()` into the directory
where this README file is located. Then you need to run the following two
lines of R code:

```{r}
	library('ProjectTemplate')
	load.project()
```

> After you enter the second line of code, you'll see a series of automated
messages as ProjectTemplate goes about doing its work. This work involves:
* Reading in the global configuration file contained in `config`.
* Loading any R packages you listed in the configuration file.
* Reading in any datasets stored in `data` or `cache`.
* Preprocessing your data using the files in the `munge` directory.  
Once that's done, you can execute any code you'd like. For every analysis
you create, we'd recommend putting a separate file in the `src` directory.
If the files start with the two lines mentioned above:

```{r}
	library('ProjectTemplate')
	load.project()
```
> You'll have access to all of your data, already fully preprocessed, and
all of the libraries you want to use.
For more details about ProjectTemplate, see http://projecttemplate.net
