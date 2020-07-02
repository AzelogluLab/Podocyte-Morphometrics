dir = "G:\\08-24-17_msPodo_Caspase_Imaging_1\\08-24-2018_MsPodo_TKI_glassplate\\08-24-2018_MsPodo_TKI_glassplate_1";
output_dir = "G:\\08-24-17_msPodo_Caspase_Imaginge_1\\output_folder";
row = newArray("A", "B", "C", "D", "E", "F", "G", "H");  
col = newArray("04", "05", "06");
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
			DAPI = well+"(fld "+fld+" wv DAPI - DAPI).tif";
			//Cy3 =  well+"(fld "+fld+" wv Cy3 - Cy3).tif";
			Cy5 =  well+"(fld "+fld+" wv Cy5 - Cy5).tif";
			FITC =  well+"(fld "+fld+" wv FITC - FITC).tif";
						
			open(dir+"\\"+DAPI); 
			run("8-bit"); 
			saveAs("Tiff", output_dir+"\\DAPI-"+DAPI); 
			close();
			
			open(dir+"\\"+FITC); //Open Actn4 channel, save orig copy and enhance other copy
			run("Duplicate...", " "); 
			run("8-bit"); 
			saveAs("Tiff", output_dir+"\\CCASP-"+FITC); 
			close();
			run("Subtract Background...", "rolling=50");
			run("Enhance Contrast...", "saturated=0.3 normalize equalize"); 

			
			open(dir+"\\"+Cy5); //open actin channel, save orig copy and enhance other copy 
			run("Duplicate...", " "); 
			run("8-bit"); 
			saveAs("Tiff", output_dir+"\\ACTIN-"+Cy5); 
			close();
			run("Subtract Background...", "rolling=50");
			run("Enhance Contrast...", "saturated=0.3 normalize equalize"); 

			run("Images to Stack", "name=Stack title=[] use"); //create combined cyto channel 
			run("Z Project...", "projection=Median");
			run("Median...", "radius=2");
			run("8-bit");
			saveAs("Tiff", output_dir+"\\CYTO-"+FITC);
			close();
			close();

			
		}
	}
}
