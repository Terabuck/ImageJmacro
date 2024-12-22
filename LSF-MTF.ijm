path="C:\\Documents and Settings\\user\\Desktop\\MTFtest";
outpath="C:\\Documents and Settings\\user\\Desktop\\output_MTF";
extname =".dcm"

dbg_val=0;

Dialog.create("\tLSF - GetPaths\t");
Dialog.addString("Input Path:", path, 100);
Dialog.addString("Output Path:", outpath, 100);
Dialog.addString("Filt Name:", extname, 30);
Dialog.addNumber("Debug:", dbg_val,0,4,"");

Dialog.show();
path=Dialog.getString();
outpath=Dialog.getString();
extname=Dialog.getString();
dbg_val=Dialog.getNumber();

if(!endsWith(path,"\\"))
{
	path = path+"\\";
}
if(!endsWith(outpath,"\\"))
{
	outpath = outpath+"\\";
}

list = getFileList(path);
//AEC variables
mm_array = newArray(list.length);
angle_array = newArray(list.length);
Exp_mode_array = newArray(list.length);
Anode_array = newArray(list.length);
Filter_array = newArray(list.length);
kV_array = newArray(list.length);
mAs_array = newArray(list.length);
uAs_array = newArray(list.length);
exp_time_array = newArray(list.length);
mA_array = newArray(list.length);
MGD_array = newArray(list.length);
ESAK_array = newArray(list.length);
HVL_array = newArray(list.length);
ImageType_array = newArray(list.length);
mAs_tot_array = newArray(list.length);

//LSF variable
wire_angle_array = newArray(list.length);
lsf_array = newArray(list.length);
pxdim=0.085;
pxspacing=0.0814;
scale = 2*sqrt(2*log(2));

ROIsize = 200; //dimension of the ROI for the wire detection
Wire_detection_ROI = 200; //ROI dimension for tungsten wire detection 

for (k=0; k<list.length; k++) //Main loop: calculation performed on every image present in the "Input path" directory
{
	if(endsWith(list[k], extname))
	{
		filename = list[k];
		open(path+filename); //open DICOM images
		nThick = nSlices;
		Stat_Slice = newArray(nThick);
		Slice_index = newArray(nThick);
		nWidth = getWidth(); //width of the images
		nHeight = getHeight(); //height of the images
		nbit = bitDepth(); //bit of images
		
		//Getting DICOM info
		angle = getNumericTag("1271,1078");
		Exp_mode = getTag("0018,7060");
		mm = getNumericTag("0018,11A0");
		Anode = getTag("0018,1191");
		Filter = getTag("0018,7050");
		kV = getNumericTag("0018,0060");
		uAs = getNumericTag("0018,1153");
		mAs = getNumericTag("0018,1152");
		mAs_tot = getNumericTag("0018,9332");
		exp_time = getNumericTag("0018,1150");
		mA = getNumericTag("0018,1151");
		Pixel_Spacing = getNumericTag("0028,0030");
		MGD = getNumericTag("0040,0316");
		ESAK = getNumericTag("0040,8302"); 
		HVL = getNumericTag("0040,0314");
		SoftwareVersion = getTag("0018,1020");
		GiottoSerialNumber = getNumericTag("0018,1000");
		MachineModel = getTag("0008,1090");
		DetectorID = getTag("0018,700A");
		LastCalibDate = getNumericTag("0018,700C");
		BitStored = getNumericTag("0028,0101");
		StudyDescription = getTag("0008,1030");
		SOP_class_UID = getTag("0008,0016");
		Projection = getNumericTag("1271,1082");
		Grid = getTag("0018,1166");
		SID = getNumericTag("0018,1110");
		SOD = getNumericTag("1271,1055");
		SAD = getNumericTag("1271,1070");
		FSP_X = getNumericTag("1271,1075");
		FSP_Y = getNumericTag("1271,1074");
		Distance_Source_to_patient = getNumericTag("0018,1111");
		
		if(SOP_class_UID == " 1.2.840.10008.5.1.4.1.1.13.1.3")
		{
			ImageType =3;//RECO
		}
		else if(SOP_class_UID == " 1.2.840.10008.5.1.4.1.1.1.2.1 ")
		{
			if (Projection == 1)
			{
				ImageType = 2;//PROJECTION
			}
			else
			{
				ImageType = 1;//MAMMO
			}
		}
		else
		{
			ImageType = 0;//UNKNOW
		}
		
		mm_array[k] = mm;
		angle_array[k] = angle;
		Exp_mode_array[k] = Exp_mode;
		Anode_array[k] = Anode;
		Filter_array[k] = Filter;
		kV_array[k] = kV;
		mAs_array[k] = mAs;
		uAs_array[k] = uAs;
		mAs_tot_array[k] = mAs_tot;
		exp_time_array[k] = exp_time;
		mA_array[k] = mA;
		MGD_array[k] = MGD;
		ESAK_array[k] = ESAK;
		HVL_array[k] = HVL;
		ImageType_array[k] = ImageType;
		
		if(ImageType_array[k] != 2) //not for Projections
		{
			//FOCUSED SLICE DETECTION
			//Detection of the minimum of the standard deviation in the ROI
			selectWindow(filename);
			if (ImageType_array[k] == 3) //RECO
			{
				nWireOffset_Y = 25/0.085;
				X_Start = nWidth/2-ROIsize/2;
				Y_Start = nHeight/2-ROIsize/2-nWireOffset_Y;

				for (s=1; s<=nThick; s++)
				{
					setSlice(s);
					makeRectangle(X_Start, Y_Start, ROIsize, ROIsize);
					roiManager("Add");
					roiManager("Select", roiManager("count")-1);
					roiManager("Rename", "WireDetectionROI");
					getRawStatistics(area, mean, min, max, std);
					Stat_Slice[s-1]= max;
					Slice_index[s-1] = s;
				}
				//Research of the maximum for max value for focused slice
				max_stat_slice = 0;
				max_slice_index = 0;
				for (i=0; i<Stat_Slice.length; i++)
				{
					if(Stat_Slice[i] > max_stat_slice)
					{
						max_stat_slice = Stat_Slice[i];
						max_slice_index = i;
					}
				}
				focused_slice = max_slice_index+1;
			}
			else if(ImageType_array[k] == 1) //MAMMO
				{
					nWireOffset_X = 60/0.085;
					nWireOffset_Y = 10/0.085;
					X_Start = nWidth-ROIsize/2-nWireOffset_X;
					Y_Start = nHeight/2-ROIsize/2-nWireOffset_Y;
				}
			//LEFT-RIGHT DIRECTION
			//Detection of the wire
			selectWindow(filename);
			if (ImageType_array[k] == 3) setSlice(focused_slice); //Focused Slice
			run("Enhance Contrast", "saturated=0.35");
			makeRectangle(X_Start, Y_Start, ROIsize, ROIsize);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "ROI_LSF");
			run("Copy");
			newImage("copy_image", "16-bit Black", Wire_detection_ROI, Wire_detection_ROI, 1);
			run("Paste");
			if (ImageType_array[k] == 1) run("Invert"); //For MAMMO images
			run("Enhance Contrast", "saturated=0.35");
			run("Select All");
			getRawStatistics(area, mean, min, max, std);
			Mean_Wire_detection_ROI = mean;
			//run("Subtract...", "value=Mean_Wire_detection_ROI");
			run("Fast Filters", "  filter=[background from median] x=10 y=0 preprocessing=smooth ");
			run("Enhance Contrast", "saturated=0.35");
			
			//Detection of the peak and coordinates calculation
			x = newArray(ROIsize);
			Peak = newArray(ROIsize);
			y = newArray(ROIsize); 
			for(i = 0; i < ROIsize; i++)
			{
				x[i] = ROIsize-i;
				makeLine(i,0,i,ROIsize);
				ydistr = getProfile();
				for(j=0; j < ydistr.length; j++)
				{
					if(ydistr[j] > Peak[i])
					{
						Peak[i] = ydistr[j];
						y[i] = j;
					}
				}
				if (dbg_val ==1) print(x[i]);
			}
			if (dbg_val ==1)
			{
				selectWindow("Log");
				saveAs("Text", outpath+"Coordinates wire Left-Right direction");
				selectWindow("Log");
				run("Close");
			}
			Fit.doFit("straight line", x, y);
			m = Fit.p(1);
			q = Fit.p(0);
			R_square_lin = Fit.rSquared;
			angle_wire = -atan(Fit.p(1));
			angle_wire_degree = atan(Fit.p(1))*180/PI;
			selectWindow("copy_image");
			run("Close");
			
			
			//PRE-SAMPLING LSF
			run("Line Width...", "line=1");
			startpt = newArray(50, 80, 110, 140); //Starting point for the 4 LSF profiles
			lsfl = 32; //Half of LSF length
			lsf = newArray(2*lsfl+1);
			zmax = 4; //number of the pre-sampling LSF
			LSF_OC_stat_a = newArray(zmax);
			LSF_OC_stat_b = newArray(zmax);
			LSF_OC_stat_c = newArray(zmax);
			LSF_OC_stat_d = newArray(zmax);
			LSF_OC_stat_FHWM = newArray(zmax);
			mtf = newArray(4096);
			scale = 2*sqrt(2*log(2)); //scale factor for the FWHM calculation
			xoff = newArray(1/tan(abs(angle_wire)));
			xOC = newArray(lsf.length*xoff.length);
			lsfOC = newArray(lsf.length*xoff.length);
			lsf_chest_nipple_sum = newArray(129);
			lsf_chest_nipple_media = newArray(129);
			wire_angle_array[k] = angle_wire_degree;
			
			//Calculation of the Pre-sampling LSF
			selectWindow(filename);
			if (ImageType_array[k] == 2) setSlice(focused_slice); //Focused Slice
			makeRectangle(X_Start, Y_Start, ROIsize, ROIsize);
			run("Enhance Contrast", "saturated=0.35");
			roiManager("Add");
			roiManager("Select", roiManager("count")-1); 
			roiManager("Rename", "Pre_samp_LSF_ROI");
			run("Copy");
			newImage("Pre_sampling_LSF_ROI", "16-bit Black", Wire_detection_ROI, Wire_detection_ROI, 1);
			run("Paste");
			if (ImageType_array[k] == 1) run("Invert");
			run("Enhance Contrast", "saturated=0.35");
			selectWindow("Pre_sampling_LSF_ROI");
				
			//Alignment of profiles
			for(z=0; z<zmax;z++) 
			{
				xstart = startpt[z]*m+q;
				for(l=0; l<xoff.length;l++) 
				{
					xoff[l] = l*tan(angle_wire);
					makeLine(startpt[z]+l,floor(xstart)-lsfl,startpt[z]+l,floor(xstart)+lsfl);
					roiManager("Add");
					roiManager("Select", roiManager("count")-1); 
					roiManager("Rename", "Pre_samp_LSF");
					lsf = getProfile();
					for (b=0; b<lsf.length; b++) 
					{
						xOC[b+l*lsf.length]= (b-xoff[l])*0.085;
						lsfOC[b+l*lsf.length]=lsf[b];
					}
				}
				Delta_x = abs(m)*0.085;
				Nyquist_Presampled_Frequency = 1/(2*Delta_x);
				profiles_number = xoff.length;

				//Gaussian fit of the raw pre-sampling lsf
				Fit.doFit("gaussian", xOC,lsfOC);
				a = Fit.p(0);
				b = Fit.p(1);
				c = Fit.p(2);
				d = Fit.p(3);

				if(dbg_val == 1 && ImageType_array[k] == 1)
				{
					print("a\t"+a);
					print("b\t"+b);
					print("mu\t"+c);
					print("sigma\t"+d);
					for(j=0; j<lsfOC.length; j++)
					{
						print(d2s(xOC[j],2)+"\t"+lsfOC[j]);

					}
					selectWindow("Log");
					saveAs("Text", outpath+"Pre_sampling_LSF_raw_MAMMO"+z+1+"lsf"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(dbg_val == 1 && ImageType_array[k] == 3)
				{
					print("a\t"+a);
					print("b\t"+b);
					print("mu\t"+c);
					print("sigma\t"+d);
					for(j=0; j<lsfOC.length; j++)
					{
						print(d2s(xOC[j],2)+"\t"+lsfOC[j]);

					}
					selectWindow("Log");
					saveAs("Text", outpath+"Pre_sampling_LSF_raw_COMBO_MAMMO"+z+1+"lsf"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}				
				LSF_OC_stat_a[z] = a;
				LSF_OC_stat_b[z] = b;
				LSF_OC_stat_c[z] = c;
				LSF_OC_stat_d[z] = abs(d);
				LSF_OC_stat_FHWM[z] = abs(scale*d);					//vector containing the FWHM value of the pre-sampling lsf
				
				//Ground subtraction
				Array.getStatistics(lsfOC, min, max, mean, std);
				lsf_mean = mean;
				for(j=0; j<lsfOC.length; j++)
				{
					lsfOC[j] = lsfOC[j]-lsf_mean;
				}
				
				//normalization
				Array.getStatistics(lsfOC, min, max, mean, std);
				lsf_max = max;
				
				for(j=0; j<lsfOC.length; j++)
				{
					lsfOC[j] = lsfOC[j]/max;
				}
				
				//Calculation of the Gaussian fit
				Fit.doFit("gaussian (no offset)", xOC,lsfOC);
				a =Fit.p(0);
				b =Fit.p(1);
				c =Fit.p(2);

				//0 value for the pixel under 15 % of the peak and gaussian fit value for the tails
				for(j=0; j<lsfOC.length; j++)
				{
					if(lsfOC[j] <= 0.15)
					{
						lsfOC[j] = 0; //Put at 0 the ground value
					}
					if(lsfOC[j] > 0.15 && lsfOC[j] < 0.25)
					{
						lsfOC[j] = a*exp(-(xOC[j]-b)*(xOC[j]-b)/(2*c*c)); //Gaussian fit of the lsf tails;
					}
				}
				
				//sort of the array
				xOC_sort = Array.copy(xOC);
				Array.sort(xOC_sort);
				x_OC_rank = Array.rankPositions(xOC);
				lsfOC_sort = newArray(lsfOC.length);
				for(j=0; j<lsfOC_sort.length; j++)
				{
					index = x_OC_rank[j];
					lsfOC_sort[j] = lsfOC[index];
				}
				
				//Re-sampling of the lsf to obtain a vector length equal to a power of 2 for fft calculation
				bFound = 0;
				nbit = 16;
				while(!bFound)
				{
					if(xOC_sort.length >= (1 << nbit))
					{
						slice_number = (1 << (nbit - 1));
						bFound = 1;
					}
					nbit --;
				}
				
				x_slice = Array.slice(xOC_sort, xOC_sort.length/2- slice_number, xOC_sort.length/2+ slice_number); //re-sampling with power of two for MTF calcultaion
				lsf_slice = Array.slice(lsfOC_sort, lsfOC_sort.length/2- slice_number, lsfOC_sort.length/2+ slice_number);
				lsf_buffer = Array.copy(lsf_slice);
				lsf_to_print = Array.copy(lsf_buffer);
				
				maxindex = b/Delta_x;
				
				//Alignment of the final lsf for printing

				lsf_buffer = Array.slice(lsf_buffer,maxindex-128,maxindex+128);
				lsf_to_print = Array.slice(lsf_to_print,maxindex-128,maxindex+128);
				for(i=0; i<lsf_buffer.length; i++)
				{
					lsf_to_print[i] += lsf_buffer[i];
					if(dbg_val == 1) print(lsf_buffer[i]);
				}
				if(ImageType_array[k] == 1 && dbg_val == 1) //MAMMO
				{
					print("a\t"+a);
					print("mu\t"+b);
					print("sigma\t"+c);
					selectWindow("Log");
					saveAs("Text", outpath+"Pre_sampling_LSF_normalized_buffer_MAMMO"+z+1+"lsf"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(ImageType_array[k] == 3 && dbg_val == 1) //RECO
				{
					print("a\t"+a);
					print("mu\t"+b);
					print("sigma\t"+c);
					selectWindow("Log");
					saveAs("Text", outpath+"Pre_sampling_LSF_normalized_buffer_RECO"+z+1+"lsf"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				
				//LSF for FFT
				//Creation of the complex array with length 2n for fft algorithm
				//real part = even element; imm part = odd element;
				signal = Array.copy(lsf_slice);
				signal_for_fft = newArray((signal.length)*2);
				for(i=0; i<(signal_for_fft.length)/2; i++)
				{
					signal_for_fft[2*i] = signal[i];
					signal_for_fft[2*i+1] = 0;
				}
				array_dim = lsf_slice.length;
				isign = 1; //For direct fft
				
				//Fourier Transform of the lsf signal
				four1(signal_for_fft,array_dim, 1);

				modulo = newArray((signal_for_fft.length)/2);
				modulo_norm = newArray((signal_for_fft.length)/2);
				mtf = Array.trim(mtf,(signal_for_fft.length)/2);
				for(i=0; i<modulo.length; i++)
				{
					modulo[i] = sqrt(signal_for_fft[2*i]*signal_for_fft[2*i]+signal_for_fft[2*i+1]*signal_for_fft[2*i+1]); //Calculation of the module of the output complex array;
				}
				Array.getStatistics(modulo, min, max, mean, std);
				norm_factor = max;
				for(i=0; i<modulo_norm.length; i++)
				{
					modulo_norm[i] = (modulo[i])/ norm_factor; //normalization
					mtf[i] += modulo_norm[i];
				}
				for (i =0; i< modulo_norm.length; i++)
				{
					if(dbg_val == 1) print(modulo_norm[i]);
				}
				if(dbg_val == 1 && ImageType_array[k] == 1)
				{						
					selectWindow("Log");
					saveAs("Text", outpath+"MTF_MAMMO"+(z+1)+"_"+mm+"mm"+d2s(mAs,0)+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(dbg_val == 1 && ImageType_array[k] == 3)
				{						
					selectWindow("Log");
					saveAs("Text", outpath+"MTF_COMBO_MAMMO"+(z+1)+"_"+mm+"mm"+d2s(mAs,0)+"mAs");
					selectWindow("Log");
					run("Close");
				}					
			} //End of the loop for the lsf acquisition
			
			x_to_print = newArray(lsf_buffer.length);
			for(i=0; i<lsf_buffer.length; i++)
			{
				x_to_print[i] = i*Delta_x;
				lsf_to_print[i] = lsf_to_print[i]/zmax;
			}
			Array.getStatistics(lsf_to_print, min, max, mean, std);
			max_lsf_to_print = max;
			
			for(i=0; i<lsf_buffer.length; i++)
			{
				lsf_to_print[i] = lsf_to_print[i]/max_lsf_to_print;
			}
			Fit.doFit("gaussian (no offset)", x_to_print,lsf_to_print);
			a_lsf_to_print =Fit.p(0);
			b_lsf_to_print =Fit.p(1);
			c_lsf_to_print =Fit.p(2);
			
			if(ImageType_array[k] == 1 && dbg_val == 1) //MAMMO
			{
				for(i=0; i<lsf_buffer.length; i++)
				{
					print(d2s(x_to_print[i],3)+"\t"+d2s(lsf_to_print[i],2)+"\t"+d2s(a_lsf_to_print*exp(-((x_to_print[i]-b_lsf_to_print)*(x_to_print[i]-b_lsf_to_print))/(2*c_lsf_to_print*c_lsf_to_print)),2));
				}
				print("a\t"+d2s(a_lsf_to_print,2));
				print("b\t"+d2s(b_lsf_to_print,2));
				print("c\t"+d2s(c_lsf_to_print,2));
				selectWindow("Log");
				saveAs("Text", outpath+"LSF_MAMMO"+mm+"mm"+d2s(mAs,0)+"mAs");
				selectWindow("Log");
				run("Close");
			}
			if(ImageType_array[k] == 3 && dbg_val == 1) //RECO
			{
				for(i=0; i<lsf_buffer.length; i++)
				{
					lsf_to_print[i] = lsf_to_print[i]/zmax;
					print(d2s(x_to_print[i],3)+"\t"+d2s(lsf_to_print[i],2)+"\t"+d2s(a_lsf_to_print*exp(-((x_to_print[i]-b_lsf_to_print)*(x_to_print[i]-b_lsf_to_print))/(2*c_lsf_to_print*c_lsf_to_print)),2));
				}
				print("a\t"+d2s(a_lsf_to_print,2));
				print("b\t"+d2s(b_lsf_to_print,2));
				print("c\t"+d2s(c_lsf_to_print,2));
				selectWindow("Log");
				saveAs("Text", outpath+"LSF_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot");
				selectWindow("Log");
				run("Close");
			}
			mtf = Array.trim(mtf, (mtf.length)/2);
			x_mtf = newArray(mtf.length);
			frequency_factor = 	Nyquist_Presampled_Frequency/((mtf.length)-1); //over sampling Nyquist frequency;
			fN_real = 1/(2*0.085); //real Nyquist frequency;
			for(i=0; i<mtf.length; i++)
			{
				x_mtf[i] = i*frequency_factor;
				mtf[i] = (mtf[i])/zmax;
				
			}
			slice_factor = fN_real/frequency_factor;
			mtf_to_print = Array.trim(mtf,slice_factor+2);

			print("mm^-1\tMTF(%)");
			for(i=0; i<mtf_to_print.length; i++)
			{
				print(d2s(x_mtf[i],2)+"\t"+d2s(mtf[i],2));
			}
			if(ImageType_array[k] == 1)//MAMMO
			{					
				selectWindow("Log");
				saveAs("Text", outpath+"MTF_MAMMO"+mm+"mm"+d2s(mAs,0)+"mAs");
				selectWindow("Log");
				run("Close");
			}
			if(ImageType_array[k] == 3)//RECO
			{					
				selectWindow("Log");
				saveAs("Text", outpath+"MTF_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot");
				selectWindow("Log");
				run("Close");
			}				

			//Statistic of the gaussian fit parameters of Pre-sampled lsf raw
			Array.getStatistics(LSF_OC_stat_a, min, max, mean, stdDev);
			mean_a_LSF = mean;
			Array.getStatistics(LSF_OC_stat_b, min, max, mean, stdDev);
			mean_b_LSF = mean;
			Array.getStatistics(LSF_OC_stat_c, min, max, mean, stdDev);
			mean_c_LSF = mean;
			Array.getStatistics(LSF_OC_stat_d, min, max, mean, stdDev);
			mean_sigma_LSF = mean;
			dev_sigma = stdDev;
			Array.getStatistics(LSF_OC_stat_FHWM, min, max, mean, stdDev);
			mean_FWHM_LSF = mean;
			dev_FWHM = stdDev;
			lsf_array[k] = mean_FWHM_LSF;

			print("Linear Fit parameters:\t");
			print("m:\t"+d2s(m,3));
			print("q:\t"+d2s(q,3));
			print("R^2:\t"+d2s(R_square_lin,3));
			print("angle (degree):\t"+d2s(angle_wire_degree,3));
			print("Over sampling pixel dim. (mm):\t"+d2s(Delta_x,4));
			print("Over sampling Nyquist frequency (mm^-1):\t"+d2s(Nyquist_Presampled_Frequency,2));
			print("Pixel dim. (mm):\t0.085");
			print("Nyquist frequency (mm^-1):\t"+d2s(fN_real,2));
			print("Used"+profiles_number+"profiles for every Over-sampled LSF");
			print("Gaussian fit parameters for Pre-sampled LSF");
			print("a\tb\tmu\tsigma");
			
			for(j=0; j<zmax; j++)
			{
				print(d2s(LSF_OC_stat_a[j],2)+"\t"+d2s(LSF_OC_stat_b[j],2)+"\t"+d2s(LSF_OC_stat_c[j],2)+"\t"+d2s(LSF_OC_stat_d[j],3));
			}
			print("Mean values\na\tb\tmu\tsigma");
			print(d2s(mean_a_LSF,2)+"\t"+d2s(mean_b_LSF,2)+"\t"+d2s(mean_c_LSF,2)+"\t"+d2s(mean_sigma_LSF,3));
			print("FWHM:\t"+d2s(mean_FWHM_LSF,3)+"\t +/-"+d2s(dev_FWHM,3));
			
			if(ImageType_array[k] == 1)
			{
				selectWindow("Log");
				saveAs("Text", outpath+"Pre_sampling_LSF_Fit_Results_MAMMO"+mm+"mm"+mAs+"mAs");
				selectWindow("Log");
				run("Close");
			}
			if(ImageType_array[k] == 3)
			{
				selectWindow("Log");
				saveAs("Text", outpath+"Pre_sampling_LSF_Fit_Results_RECO"+mm+"mm"+mAs_tot+"mAstot");
				selectWindow("Log");
				run("Close");
			}
		}	//End of the if image is not a projection
	}
if(dbg_val ==1 && ImageType_array[k] == 1 ) roiManager("Save", outpath+"RoiSet_MAMMO"+mm+"mm"+mAs+"mAs"+".zip");
if(dbg_val ==1 && ImageType_array[k] == 3) roiManager("Save", outpath+"RoiSet_RECO"+mm+"mm"+mAs_tot+"mAstot"+".zip");

run("Close All");
roiManager("Reset");
selectWindow("ROI Manager");
run("Close");

}	//End of the Main loop performed on images present in Input path directory

//Printing of the LSF-MTF Test Results
print("LSF-MTF Test");
print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
print("Machine Model:\t"+MachineModel);
print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
print("Detector ID:\t"+DetectorID);
print("SoftwareVersion:\t"+SoftwareVersion);
print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
print("\nmmPMMA\tFILTER\tkV\tmAs\tangle (degree)\tFHWM (mm)");
for(j = 0; j < list.length; j++)
{
	if (ImageType_array[j] == 1)//MAMMO
	{
		print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(wire_angle_array[j],3)+"\t"+d2s(lsf_array[j],3));
	}
}
selectWindow("Log");
saveAs("Text", outpath+"LSF_Result_MAMMO.txt");
selectWindow("Log");
run("Close");

//RECO IMAGES RESULTS
print("LSF-MTF Test");
print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
print("Machine Model:\t"+MachineModel);
print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
print("Detector ID:\t"+DetectorID);
print("SoftwareVersion:\t"+SoftwareVersion);
print("Last Calibration Date:\t"+d2s(LastCalibDate,0));

print("\nmmPMMA\tFILTER\tkV\tmAs\tangle (degree)\tFHWM (mm)");
for(j = 0; j < list.length; j++)
{
	if (ImageType_array[j] == 3)//RECO
	{
		print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(wire_angle_array[j],3)+"\t"+d2s(lsf_array[j],3));
	}
}
selectWindow("Log");
saveAs("Text", outpath+"LSF_Result_RECO.txt");
selectWindow("Log");
run("Close");

//FFT function from Numerical Recipes 

function four1(data, n, isign) 
{
	
	/*Replaces data[0..2*n-1] by its discrete Fourier transform, if isign is input as 1; or replaces
	data[0..2*n-1] by n times its inverse discrete Fourier transform, if isign is input as -1. 
	data is a complex array of length n stored as a real array of length 2*n. 
	n must be an integer power of 2.*/
	//Int nn,mmax,m,j,istep,i;
	//Doub wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
	//if (n<2 || n&(n-1)) print("n must be power of 2 in four1");
	
	nn = n*2;
	j = 1;
	
	for (i=1;i<nn;i+=2) //This is the bit-reversal section of the routine.
	{
		if (j > i) 
		{ 
			temp1 = data[j-1];
			data[j-1] = data[i-1];
			data[i-1] = temp1;
			
			temp2 = data[j];
			data[j] = data[i];
			data[i] = temp2;
		}
		m = n;
		while (m >= 2 && j > m) 
		{
			j -= m;
			m = m/2;
		}
		j += m;
	}
	// Here begins the Danielson-Lanczos section of the routine.
	mmax=2;
	while (nn > mmax) 								//Outer loop executed log2 n times.
	{ 
		istep = mmax *2;
		theta = isign*(2*PI/mmax); 		//Initialize the trigonometric recurrence.
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m +=2) 				//Here are the two nested inner loops.
		{ 
			for (i = m; i<=nn; i+=istep) 
			{
				j = i+mmax; 						//This is the Danielson-Lanczos formula:
				tempr = wr*data[j-1]-wi*data[j]; 
				tempi = wr*data[j]+wi*data[j-1];
				data[j-1] = data[i-1]-tempr;
				data[j] = data[i]-tempi;
				data[i-1] += tempr;
				data[i] += tempi;
			}
			wtemp = wr;
			wr = wtemp*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax=istep;
	}
}

//Get DICOM INFO
function getNumericTag(tag) 
{
    value = getTag(tag);
    if (value=="") return NaN;
    index3 = indexOf(value, "\\");
    if (index3>0) value = substring(value, 0, index3);
    value = 0 + value; // convert to number
    return value;
 }

function getTag(tag) 
{
      info = getImageInfo();
      index1 = indexOf(info, tag);
      if (index1==-1) return "";
      index1 = indexOf(info, ":", index1);
      if (index1==-1) return "";
      index2 = indexOf(info, "\n", index1);
      value = substring(info, index1+1, index2);
      return value;
}
