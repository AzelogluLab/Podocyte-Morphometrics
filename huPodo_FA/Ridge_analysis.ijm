dir = "F:\\09-04-2018_Dose_Response_FA\\08-24-2018_MsPodo_TKI_glassplate\\08-24-2018_MsPodo_TKI_glassplate_2";
output_dir = "F:\\09-04-2018_Dose_Response_FA\\output_ridge";
row = newArray("A", "B", "C", "D", "E", "F", "G", "H");  
col = newArray("10", "11", "12");
setBatchMode(true); 

for (r=0; r<lengthOf(row); r++){
	for (c=0; c<lengthOf(col); c++){
		for (f=1; f<=36; f++){
			well = row[r]+" - "+col[c];
			if (f<10){
				fld = "0"+f;
			} else {
				fld = f; 
			}
			ACTIN =  well+"(fld "+fld+" wv Cy5 - Cy5).tif";
						
			open(dir+"\\"+ACTIN);
			run("Subtract Background...", "rolling=50");
			setMinAndMax(86, 3902);

			run("Apply LUT"); 
			run("8-bit");
			run("Ridge Detection", "line_width=1 high_contrast=250 low_contrast=80 estimate_width make_binary method_for_overlap_resolution=NONE sigma=0.60 lower_threshold=5 upper_threshold=18 minimum_line_length=0 maximum=0");
			selectWindow(ACTIN + " Detected segments"); 
			run("Invert LUT"); 
			rename("RIDGES-"+ACTIN); 
			saveAs("Tiff", output_dir+"\\"+getTitle());
			close();
			close();

				
		}
	}
}
