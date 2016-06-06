setPaths <- function(dat_scenario, iter, Fdyn, Rdyn, SigmaR, 
	test_dir, iterpath, SimDir, Ncomp=NULL, lh, sens_param, sens_label,
    RecTypeDesc, OMselex, EMselex){

    simdir <- file.path(test_dir, SimDir)
    dir.create(simdir, showWarnings=FALSE)

    if(SimDir=="base"){
    	lhdir <- file.path(simdir, paste0("LH_", lh))
    	dir.create(lhdir, showWarnings=FALSE)
    }

    if(SimDir=="explore_comp"){
    	compdir <- file.path(simdir, paste0("Nyrs_comp_", Ncomp))
    	dir.create(compdir, showWarnings=FALSE)

    	lhdir <- file.path(compdir, paste0("LH_", lh))
    	dir.create(lhdir, showWarnings=FALSE)
    }

    if(SimDir=="sensitivity"){
        paramdir <- file.path(simdir, paste0(sens_param, "_", sens_label))
        dir.create(paramdir, showWarnings=FALSE)

        lhdir <- file.path(paramdir, paste0("LH_", lh))
        dir.create(lhdir, showWarnings=FALSE)
    }

    selexdir <- file.path(lhdir, paste0("OM", OMselex,"_EM", EMselex))
    dir.create(selexdir, showWarnings=FALSE)

    recdir <- file.path(selexdir, paste0("RecType_", RecTypeDesc))
    dir.create(recdir, showWarnings=FALSE)

    sigrdir <- file.path(recdir, paste0("SigR_", SigmaR))
    dir.create(sigrdir, showWarnings=FALSE)

	Fdir <- file.path(sigrdir, paste0("Fdyn_", Fdyn))
	dir.create(Fdir, showWarnings=FALSE)

	Rdir <- file.path(Fdir, paste0("Rdyn_", Rdyn))
	dir.create(Rdir, showWarnings=FALSE)

	moddir <- file.path(Rdir, dat_scenario)
	dir.create(moddir, showWarnings=FALSE)

	if(iterpath==TRUE){
	  iterdir <- file.path(moddir, iter)
	  dir.create(iterdir, showWarnings=FALSE)
	}

	if(iterpath==TRUE) return(iterdir)
    if(iterpath==FALSE) return(moddir)

}