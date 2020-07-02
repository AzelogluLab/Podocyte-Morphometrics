dir = "G:\\HCA_trial_9_5_2018";
output_dir = "G:\\HCA_trial_9_5_2018\\output";
row = newArray("A", "B", "C", "D", "E", "F", "G", "H");  
col = newArray("01", "02", "03");
setBatchMode(true); 

for (r=0; r<lengthOf(row); r++){
	for (c=0; c<lengthOf(col); c++){
		for (f=1; f<=7; f++){
			well = row[r]+" - "+col[c];
			if (f<10){
				fld = "0"+f;
			} else {
				fld = f; 
			}
			DAPI = well+"(fld "+fld+" wv DAPI - DAPI).tif";
			Cy3 =  well+"(fld "+fld+" wv Cy3 - Cy3).tif";
			Cy5 =  well+"(fld "+fld+" wv Cy5 - Cy5).tif";
			FITC =  well+"(fld "+fld+" wv FITC - FITC).tif";
						
			open(dir+"\\"+Cy3);
			run("Subtract Background...", "rolling=50");
			run("CLAHE ", "blocksize=19 histogram=256 maximum=3");
			run("Enhance Contrast...", "saturated=0.35");	
			mm = getTitle();
			run("LoG 3D", "sigmax=2 sigmay=2");
			selectWindow("LoG of "+mm);
			run("8-bit");
			run("Make Binary");
			run("Invert LUT"); 
			rename("BPAX-"+Cy3);
			saveAs("Tiff", output_dir+"\\"+getTitle());
			close();  
			close();
			open(dir+"\\"+DAPI); 
			run("8-bit"); 
			saveAs("Tiff", output_dir+"\\DAPI-"+DAPI); 
			close();
			
			open(dir+"\\"+FITC); //Open Actn4 channel, save orig copy and enhance other copy
			run("Duplicate...", " "); 
			run("8-bit"); 
			saveAs("Tiff", output_dir+"\\ACTN4-"+FITC); 
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
