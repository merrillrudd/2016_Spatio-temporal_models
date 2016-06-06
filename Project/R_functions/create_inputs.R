create_inputs <- function(param, val, lh_dat){
	
		dat_input <- lh_dat
		dat_input$log_sigma_C <- log(lh_dat$SigmaC)
		dat_input$log_sigma_I <- log(lh_dat$SigmaI)
		dat_input$log_CV_L <- log(lh_dat$CVlen)
		if(param!=FALSE) dat_input[[param]] <- val

	return(dat_input)
}