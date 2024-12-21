//ImageJ Macro written by Dott. Turco Alessandro IMS Italy (BO) June-August 2013
//Macro modified by Dott. Turco Alessandro April 2014 
//MACRO VERSION 1.2.1

requires("1.47a");

outpath = "D:\\";
extname =".dcm"
pxdim = 0.085; //pixel dimension (mm)
detailpos = 60; //position of the Al object with respect to the chest edge expressed in mm
ROIsize = 500; //ROI dimension for the initial object research
Object_ROIsize = 60; //ROI dimension for the statistical analysis of CNR
context = 1;
LSF_cal = 1;
dbg_val = 0;

myDir = outpath +"AEC_Test_Results"+File.separator;
File.makeDirectory(myDir);
if (!File.exists(myDir)) exit("Unable to create directory");

Dialog.create("\tAEC Test - GetPaths\t");

Dialog.addNumber("Object Position (mm from chest side) :", detailpos,0,3,"");
Dialog.addNumber("Context (1= MAMMO; 2= TOMO; 3= COMBO):", context,0,4,"");
Dialog.addNumber("Calculate the LSF profiles (0= Not Cal; 1= Cal):", LSF_cal,0,4,"");
Dialog.addNumber("Debug:", dbg_val,0,4,"");
Dialog.show();


detailpos = Dialog.getNumber();
context = Dialog.getNumber();
/* context Legend: 
	context = 1 -->  MAMMO;
	context = 2 -->  TOMO;
	context = 3 -->  COMBO;
*/ 
LSF_cal = Dialog.getNumber();
dbg_val = Dialog.getNumber();

if(context == 1) context_dialog = "MAMMO";
if(context == 2) context_dialog = "TOMO";
if(context == 3) context_dialog = "COMBO";

dir = getDirectory("Choose the directory containing the AEC_"+context_dialog+" images");
list = getFileList(dir);

detailposition = floor(detailpos/pxdim);

run("Set Measurements...", "mean standard scientific redirect=None decimal=6");

if(dbg_val==1) print("number of files: "+d2s(list.length,0));


		//Creation of arrays used for sort of results
		
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
		
		//CNR variables
		MeanBKG_array = newArray(list.length);
		StdBKG_array = newArray(list.length);
		SNR_BKG_array = newArray(list.length);
		MeanDetail_array = newArray(list.length);
		StdDetail_array = newArray(list.length);
		SNR_Object_array = newArray(list.length);
		Contrast_array = newArray(list.length);
		CNR_array = newArray(list.length);
		result_SNR_array = newArray(list.length);
		Mean_Dev_array = newArray(5);
		SNR_Dev_array = newArray(5);
		Hor_Al_Profile_length = 500;
		BackGround_mean_array = newArray(list.length);
		BackGround_std_array = newArray(list.length);
		BackGround_SNR_array = newArray(list.length);
		Al_mean_array = newArray(list.length);
		Al_std_array = newArray(list.length);
		Al_SNR_array = newArray(list.length);
		Contrast_rec_array = newArray(list.length);
		CNR_rec_array = newArray(list.length);
		
		//Homogeneity variables
		MeanROI1 = newArray(list.length);
		MeanROI2 = newArray(list.length);
		MeanROI3 = newArray(list.length);
		MeanROI4 = newArray(list.length);
		MeanROI5 = newArray(list.length);
		SNR_ROI1 = newArray(list.length);
		SNR_ROI2 = newArray(list.length);
		SNR_ROI3 = newArray(list.length);
		SNR_ROI4 = newArray(list.length);
		SNR_ROI5 = newArray(list.length);
		Mean_Dev1 = newArray(list.length);
		Mean_Dev2 = newArray(list.length);
		Mean_Dev3 = newArray(list.length);
		Mean_Dev4 = newArray(list.length);
		Mean_Dev5 = newArray(list.length);
		SNR_Dev1 = newArray(list.length);
		SNR_Dev2 = newArray(list.length);
		SNR_Dev3 = newArray(list.length);
		SNR_Dev4 = newArray(list.length);
		SNR_Dev5 = newArray(list.length);
		Max_Mean_Dev = newArray(list.length);
		Max_SNR_Dev = newArray(list.length);
		MeanTot = newArray(list.length);
		SNRTot = newArray(list.length);
		
		//FFT variables
		NPS_ROI_size = 256;
		nROI = 16;
		
		Hor_Profile = newArray(0.5*NPS_ROI_size);
		Hor_Profile2 = newArray(0.5*NPS_ROI_size);
		Ver_Profile = newArray(0.5*NPS_ROI_size);
		Ver_Profile2 = newArray(0.5*NPS_ROI_size);
		buffer_hor = newArray(0.5*NPS_ROI_size);
		buffer_ver = newArray(0.5*NPS_ROI_size);
		buffer_rad = newArray(0.5*NPS_ROI_size);
		rad_upper = newArray(0.5*NPS_ROI_size);
		rad_lower = newArray(0.5*NPS_ROI_size);
		NPS_rad = newArray(0.5*NPS_ROI_size);
		hor_upper = newArray(0.5*NPS_ROI_size);
		hor_lower = newArray(0.5*NPS_ROI_size);
		NPS_hor = newArray(0.5*NPS_ROI_size);
		ver_upper = newArray(0.5*NPS_ROI_size);
		ver_lower = newArray(0.5*NPS_ROI_size);
		NPS_ver = newArray(0.5*NPS_ROI_size);
		grid_visibility_array = newArray(list.length);
		Grid_Peak_array = newArray(list.length);
		NPS_Ground_array = newArray(list.length);
		
		//ASF variable
		sigma_asf = newArray(list.length);
		CNR_asf_focused_slice = newArray(list.length);
		
		//LSF variable
		wire_angle_array = newArray(list.length);
		lsf_array = newArray(list.length);
		
		//Alignment Test variable
		left_right_alingment_results = newArray(list.length);
		chest_nipple_alingment_results = newArray(list.length);
		focused_slice_chest_array = newArray(list.length);
		focused_slice_nipple_array = newArray(list.length);
		focused_slice_left_array = newArray(list.length);
		focused_slice_right_array = newArray(list.length);
		SliceNumber_array = newArray(list.length);
		
		//Test Result variables
		result_SNR = "Fail";
		result_exp_time = "Fail";
		result_Homogeneity_Mean ="Fail";
		result_Homogeneity_SNR ="Fail";
		grid_visibility_result = "Fail";

for (k=0; k<list.length; k++) //Main loop performed for all the images in the "input dir" directory
{
	if(endsWith(list[k], extname))
	{
		filename = list[k];
		//for DICOM images
		if(extname==".dcm")open(dir+filename); 
		//for RAW images
		if(extname==".raw")run("Raw...", "open=["+dir+filename+"] image=[16-bit Unsigned] width=2816 height=3584 offset=0 number=1 gap=0 white little-endian");

		if(extname==".dcm")
		{
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
			
		}
		/* ImageType Legend: 
			ImageType = 0 -->  UNKNOW;
			ImageType = 1 -->  MAMMO;
			ImageType = 2 -->  PROJECTION;
			ImageType = 3 -->  TOMO RECO;
		*/ 
		if(SOP_class_UID == " 1.2.840.10008.5.1.4.1.1.13.1.3")
		{
			ImageType =3;//TOMO RECO
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
		
		nWidth = getWidth(); //width of the images
		nHeight = getHeight(); //height of the images
		nbit = bitDepth(); //bit of images
		
		// Al OBJECYT DETECTION for TOMO RECO images 
		//Detection of the focused slice algorithm 
		if (ImageType_array[k]==3) //TOMO RECO
		{
			nThick = nSlices;
			SliceNumber_array[k] = nThick;
			Stat_Slice = newArray(nThick);
			Slice_index = newArray(nThick);
			//Detection of the minimum of the standard deviation in the ROI
			selectWindow(filename);
			for (s=1; s<=nThick; s++)
			{
				setSlice(s);
				makeRectangle(nWidth-detailposition-0.5*ROIsize, 0.5*nHeight-0.5*ROIsize, ROIsize, ROIsize);
				getRawStatistics(area, mean, min, max, std, histogram);
				Stat_Slice[s-1]= std;
				Slice_index[s-1] = s;
			}
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "ROI for focused slice detection");
			//Research of the maximum for std. dev.
			max_stat_slice = 0;
			max_slice_index = 0;
			for (i=0; i<Stat_Slice.length; i++)
			{
				if(Stat_Slice[i]>max_stat_slice)
				{
					max_stat_slice = Stat_Slice[i];
					max_slice_index = i;
				}
			}
			//Al object location
			selectWindow(filename);
			focused_slice = max_slice_index+1;
			setSlice(focused_slice); //Focused Slice
			makeRectangle(nWidth-detailposition-0.5*ROIsize, 0.5*nHeight-0.5*ROIsize, ROIsize, ROIsize);
			run("Copy");
			newImage("copy_image", "16-bit Black", ROIsize, ROIsize, 1);
			run("Paste");
			run("Enhance Contrast", "saturated=0.35");
			run("Fast Filters", "link filter=median x=20 y=20 ");//median filter for the smoothing for the edge profile
			run("Line Width...", "line="+d2s(ROIsize,0));
			selectWindow("copy_image");
			
			//locate object in Y direction
			makeLine(ROIsize/2,0,ROIsize/2,ROIsize);
			if(dbg_val==1) roiManager("Add");
			
			//ydistr = vector containing the profile (double edge)
			ydistr = getProfile();
							
			//xderiv = vector containing the derivate of profile (two peaks)
			yderiv=newArray(ydistr.length-1);
			xderiv=newArray(ydistr.length-1);
			
			maxyderiv = 0;
			minyderiv = 0;
			maxindex = 0;
			minindex = 0;
			
			//locate centre in Y direction finding the maximum and minimum of the derivate
			for(i=0; i<ydistr.length-1; i++)
			{
				xderiv [i] = i;
				yderiv[i]=(ydistr[i+1]-ydistr[i])/(i+1-i);
				
				if (yderiv[i] > maxyderiv)
				{
					maxyderiv = yderiv[i];
					maxindex = i;
				}
				if (yderiv[i] < minyderiv)
				{
					minyderiv = yderiv[i];
					minindex = i;
				}
			}
			if(dbg_val == 1) 
			{
				Plot.create("Derivate in Y direction","pixel","a.u.");
				Plot.setLimits(0, ydistr.length-1, -2, 2);
				Plot.add("Deriv", xderiv, yderiv);
				Plot.show();
			}
			Object_lenght_Y = maxindex-minindex;
			Y = (maxindex+minindex)/2;
			
			if(dbg_val == 1) 
			{
				print("Y DIRECTION maxindex:\t"+d2s(maxindex,0)+"\t minindex:\t"+d2s(minindex,0));
				selectWindow("Log");
				saveAs("Text", myDir+"Y_DIRECTION_COORDINATES"+mm+"mm"+mAs+"mAs");
				selectWindow("Log");
				run("Close");
			}
			Yobject = 0.5*nHeight-0.5*ROIsize+Y;
			
			//locate object in X direction
			selectWindow("copy_image");
			run("Enhance Contrast", "saturated=0.35");
			makeLine(0,ROIsize/2,ROIsize,ROIsize/2);
			if(dbg_val == 1) roiManager("Add");
			//ydistr = vector containing the profile (double edge)
			ydistr = getProfile();
			if(dbg_val == 1) run("Plot Profile");
			//xderiv = vector containing the derivate of profile (two peaks)
			yderiv=newArray(ydistr.length-1);
			xderiv=newArray(ydistr.length-1);
			
			maxyderiv = 0;
			minyderiv = 0;
			maxindex = 0;
			minindex = 0;
			
			//locate centre in X direction finding the maximum and minimum of the derivate
			for(i=0; i<ydistr.length-1; i++)
			{
				xderiv [i] = i;
				yderiv[i]=(ydistr[i+1]-ydistr[i])/(i+1-i);
				//max value research
				if (yderiv[i] > maxyderiv)
				{
					maxyderiv = yderiv[i];
					maxindex = i;
				}
				if (yderiv[i] < minyderiv)
				{
					minyderiv = yderiv[i];
					minindex = i;
				}
			}
			if(dbg_val==1) 
			{
				Plot.create("Derivate in X direction","pixel","a.u.");
				Plot.setLimits(0, ydistr.length-1, -2, 2);
				Plot.add("Deriv", xderiv, yderiv);
				Plot.show();
			}
			Object_lenght_X = maxindex-minindex;
			X = (maxindex+minindex)/2;
			
			if(dbg_val == 1)
			{
				print("X DIRECTION maxindex:\t"+d2s(maxindex,0)+"\t minindex:\t"+d2s(minindex,0));
				selectWindow("Log");
				saveAs("Text", myDir+"X DIRECTION COORDINATES"+mm+"mm"+mAs+"mAs");
				selectWindow("Log");
				run("Close");
			}	
			
			Xobject = nWidth-detailposition-0.5*ROIsize+X;
			selectWindow("copy_image");
			run("Close");
			
			//Detection of the focused slice using tungsten wire object at chest side
			Stat_Slice_wire = newArray(nThick);
			Slice_index_wire = newArray(nThick);
			LSF_Offset_X = 35/pxdim; //Position of the tungsten wire with respect to the Al object in mm on X direction
			LSF_Offset_Y = 10/pxdim; //Position of the tungsten wire with respect to the Al object in mm on Y direction
			Wire_detection_ROI = 200; //ROI dimension for tungsten wire detection 
			
			//Detection of the maximum in the ROI
			selectWindow(filename);
			for (c=1; c<=nThick; c++)
			{
				setSlice(c);
				makeRectangle(Xobject+LSF_Offset_X-0.5*Wire_detection_ROI, Yobject-LSF_Offset_Y-0.5*Wire_detection_ROI, Wire_detection_ROI, Wire_detection_ROI);	
				getRawStatistics(area, mean, min, max, std, histogram);
				Stat_Slice_wire[c-1]= std;
				Slice_index_wire[c-1] = c;
			}
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "ROI wire chest side");
			//Research of the maximum for std. dev.
			max_stat_slice_wire = 0;
			max_slice_index_wire = 0;
			for (q=0; q<Stat_Slice_wire.length; q++)
			{
				if(Stat_Slice_wire[q]>max_stat_slice_wire)
				{
					max_stat_slice_wire = Stat_Slice_wire[q];
					max_slice_index_wire = q;
				}
			}
			focused_slice_real = max_slice_index_wire + 1;
			
			//ALIGNMENT TEST
			if(LSF_cal == 1)
			{
				//Detection of the focused slice using tungsten wire object in the nipple side
				Stat_Slice_wire_nipple = newArray(nThick);
				Slice_index_wire_nipple = newArray(nThick);
				Offset_ROI_Alignment_X = 140/pxdim; //Position of the tungsten wire with respect to the Al object in mm on X direction
				Offset_ROI_Alignment_Y = 16/pxdim;  //Position of the tungsten wire with respect to the Al object in mm on Y direction
				//Detection of the maximum in the ROI
				selectWindow(filename);
				for (c=1; c<=nThick; c++)
				{
					setSlice(c);
					makeRectangle(Xobject-Offset_ROI_Alignment_X-0.5*Wire_detection_ROI, Yobject-Offset_ROI_Alignment_Y-0.5*Wire_detection_ROI, Wire_detection_ROI, Wire_detection_ROI);
					getRawStatistics(area, mean, min, max, std, histogram);
					Stat_Slice_wire_nipple[c-1]= max;
					Slice_index_wire_nipple[c-1] = c;
				}
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI wire nipple side");
				//Research of the maximum for std. dev.
				max_stat_slice_wire_nipple = 0;
				max_slice_index_wire_nipple = 0;
				for (q=0; q<Stat_Slice_wire_nipple.length; q++)
				{
					if(Stat_Slice_wire_nipple[q]>max_stat_slice_wire_nipple)
					{
						max_stat_slice_wire_nipple = Stat_Slice_wire_nipple[q];
						max_slice_index_wire_nipple = q;
					}
				}
				focused_slice_nipple_side = max_slice_index_wire_nipple + 1;


				//Detection of the focused slice using tungsten wire object in the Right side
				Stat_Slice_wire_right = newArray(nThick-9);
				Slice_index_wire_right = newArray(nThick-9);
				Offset_ROI_Alignment_X = 55/pxdim; //Position of the tungsten wire with respect to the Al object in mm on X direction
				Offset_ROI_Alignment_Y = 120/pxdim;  //Position of the tungsten wire with respect to the Al object in mm on Y direction
				//Detection of the maximum in the ROI
				selectWindow(filename);
				for (c=10; c<=nThick; c++)
				{
					setSlice(c);
					makeRectangle(Xobject-Offset_ROI_Alignment_X-0.5*Wire_detection_ROI, Yobject-Offset_ROI_Alignment_Y-0.5*Wire_detection_ROI, Wire_detection_ROI, Wire_detection_ROI);
					getRawStatistics(area, mean, min, max, std, histogram);
					Stat_Slice_wire_right[c-10]= max;
					Slice_index_wire_right[c-10] = c;
				}
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI wire right side");
				//Detection of the maximum in the ROI
				max_stat_slice_wire_right = 0;
				max_slice_index_wire_right = 0;
				for (q=0; q<Stat_Slice_wire_right.length; q++)
				{
					if(Stat_Slice_wire_right[q]>max_stat_slice_wire_right)
					{
						max_stat_slice_wire_right = Stat_Slice_wire_right[q];
						max_slice_index_wire_right = Slice_index_wire_right[q];
					}
				}
				focused_slice_right_side = max_slice_index_wire_right;
				
				//Detection of the focused slice using tungsten wire object in the Left side
				Stat_Slice_wire_left = newArray(nThick-9);
				Slice_index_wire_left = newArray(nThick-9);
				Offset_ROI_Alignment_X = 55/pxdim; //Position of the tungsten wire with respect to the Al object in mm on X direction
				Offset_ROI_Alignment_Y = 120/pxdim;  //Position of the tungsten wire with respect to the Al object in mm on Y direction
				//Detection of the maximum in the ROI
				selectWindow(filename);
				for (c=10; c<=nThick; c++)
				{
					setSlice(c);
					makeRectangle(Xobject-Offset_ROI_Alignment_X-0.5*Wire_detection_ROI, Yobject+Offset_ROI_Alignment_Y-0.5*Wire_detection_ROI, Wire_detection_ROI, Wire_detection_ROI);
					getRawStatistics(area, mean, min, max, std, histogram);
					Stat_Slice_wire_left[c-10]= max;
					Slice_index_wire_left[c-10] = c;
				}
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI wire left");
				//Detection of the maximum in the ROI
				max_stat_slice_wire_left = 0;
				max_slice_index_wire_left = 0;
				for (q=0; q<Stat_Slice_wire_left.length; q++)
				{
					if(Stat_Slice_wire_left[q]>max_stat_slice_wire_left)
					{
						max_stat_slice_wire_left = Stat_Slice_wire_left[q];
						max_slice_index_wire_left = Slice_index_wire_left[q];
					}
				}
				focused_slice_left_side = max_slice_index_wire_left;
				
				//Alignment Test Result
				
				//Chest-Nipple Direction
				Delta_Slice_chest_nipple = abs(focused_slice_real-focused_slice_nipple_side);
				if(Delta_Slice_chest_nipple <= 1)
				{
					chest_nipple_alingment_results[k] = "Pass";
				}
				else
				{
					chest_nipple_alingment_results[k] = "Fail";
				}
				
				//Left-Right Direction
				Delta_Slice_left_right = abs(focused_slice_left_side-focused_slice_right_side);
				if(Delta_Slice_left_right <= 1)
				{
					left_right_alingment_results[k] = "Pass";
				}
				else
				{
					left_right_alingment_results[k] = "Fail";
				}
				focused_slice_chest_array[k] = focused_slice_real;
				focused_slice_nipple_array[k] = focused_slice_nipple_side;
				focused_slice_right_array[k] = focused_slice_right_side;
				focused_slice_left_array[k] = focused_slice_left_side;
			}
			/*
			//CNR Variables
			selectWindow(filename);
			setSlice(focused_slice_real);     //selection of the focused slice
			run("Enhance Contrast", "saturated=0.35");
			run("Line Width...", "line=10");
			makeLine(Xobject-Hor_Al_Profile_length/2, Yobject, Xobject+Hor_Al_Profile_length/2, Yobject); //profile of the Al object
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "CNR profile");
			Hor_Al_Profile = getProfile();
			//print("Focused Slice:"+focused_slice_real);
			for (c=0; c<=Hor_Al_Profile.length-1; c++)
			{
				print(Hor_Al_Profile[c]);
			}
			selectWindow("Log");
			saveAs("Text", myDir+"Hor_Al_Profile_slide"+focused_slice_real+"of"+nThick);
			selectWindow("Log");
			run("Close");
			*/
		}
		else //for mammo and projection images
		{
			selectWindow(filename);
			makeRectangle(nWidth-detailposition-0.5*ROIsize, 0.5*nHeight-0.5*ROIsize, ROIsize, ROIsize);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "Start ROI");
			run("Copy");
			newImage("copy_image", "16-bit Black", ROIsize, ROIsize, 1);
			run("Paste");
			run("Enhance Contrast", "saturated=0.35");
			if(context == 1) //MAMMO
			{
				run("Fast Filters", "link filter=median x=10 y=10 "); //median filter for the smoothing for the edge profile
			}
			if(context == 2 || context == 3) //TOMO or COMBO PROJ
			{
				run("Fast Filters", "filter=mean x=1 y=6 preprocessing=smooth "); //median filter for the smoothing for the edge profile
			}
			run("Line Width...", "line="+d2s(ROIsize,0));
			selectWindow("copy_image");
			
			//locate object in Y direction
			makeLine(ROIsize/2,0,ROIsize/2,ROIsize);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "Line for vertical profile");
			
			//ydistr = vector containing the profile (double edge)
			ydistr = getProfile();
			
			//xderiv = vector containing the derivate of profile (two peaks)
			yderiv=newArray(ydistr.length-1);
			xderiv=newArray(ydistr.length-1);
			
			maxyderiv = 0;
			minyderiv = 0;
			maxindex = 0;
			minindex = 0;
			
			//locate centre in Y direction finding the maximum and minimum of the derivate
			for(i=0; i<ydistr.length-1; i++)
			{
				xderiv [i] = i;
				yderiv[i]=(ydistr[i+1]-ydistr[i])/(i+1-i);
				if (yderiv[i] > maxyderiv)
				{
					maxyderiv = yderiv[i];
					maxindex = i;
				}
				if (yderiv[i] < minyderiv)
				{
					minyderiv = yderiv[i];
					minindex = i;
				}
			}
			if(dbg_val == 1) 
			{
				Plot.create("Derivate in Y direction","pixel","a.u.");
				Plot.setLimits(0, ydistr.length-1, -2, 2);
				Plot.add("Deriv", xderiv, yderiv);
				Plot.show();
			}
			Object_lenght_Y = maxindex-minindex;
			Y = (maxindex+minindex)/2;
			
			if(dbg_val == 1) 
			{
				print("Y DIRECTION maxindex:\t"+d2s(maxindex,0)+"\t minindex:\t"+d2s(minindex,0));
				selectWindow("Log");
				saveAs("Text", myDir+"Y_DIRECTION_COORDINATES"+mm+"mm"+mAs+"mAs");
				selectWindow("Log");
				run("Close");
			}
			Yobject = 0.5*nHeight-0.5*ROIsize+Y;
			
			//locate object in X direction
			selectWindow("copy_image");
			run("Enhance Contrast", "saturated=0.35");
			makeLine(0,ROIsize/2,ROIsize,ROIsize/2);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "Line for horizontal profile");
			//ydistr = vector containing the profile (double edge)
			ydistr = getProfile();

			//xderiv = vector containing the derivate of profile (two peaks)
			yderiv=newArray(ydistr.length-1);
			xderiv=newArray(ydistr.length-1);
			
			maxyderiv = 0;
			minyderiv = 0;
			maxindex = 0;
			minindex = 0;
			
			//locate centre in X direction finding the maximum and minimum of the derivate
			for(i=0; i<ydistr.length-1; i++)
			{
				xderiv [i] = i;
				yderiv[i]=(ydistr[i+1]-ydistr[i])/(i+1-i);
				//max value research
				if (yderiv[i] > maxyderiv)
				{
					maxyderiv = yderiv[i];
					maxindex = i;
				}
				if (yderiv[i] < minyderiv)
				{
					minyderiv = yderiv[i];
					minindex = i;
				}
			}
			if(dbg_val==1) 
			{
				Plot.create("Derivate in X direction","pixel","a.u.");
				Plot.setLimits(0, ydistr.length-1, -2, 2);
				Plot.add("Deriv", xderiv, yderiv);
				Plot.show();
			}
			Object_lenght_X = maxindex-minindex;
			X = (maxindex+minindex)/2;
			
			if(dbg_val==1)
			{
				print("X DIRECTION maxindex:\t"+d2s(maxindex,0)+"\t minindex:\t"+d2s(minindex,0));
				selectWindow("Log");
				saveAs("Text", myDir+"X_DIRECTION_COORDINATES"+mm+"mm"+mAs+"mAs");
				selectWindow("Log");
				run("Close");
			}	
			
			Xobject = nWidth-detailposition-0.5*ROIsize+X;
			
			selectWindow("copy_image");
			run("Close");
			selectWindow(filename);
			makeRectangle(Xobject, Yobject, 1, 1);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "Al_position");
			if(dbg_val==1) 
			{
				print("Object centre in small image");
				print("X:\t"+d2s(X,0)+"\t Y:\t"+d2s(Y,0));
				print("Object centre in real image");
				print("X real:\t"+d2s(Xobject,0)+"\t Y real:\t"+d2s(Yobject,0));
				selectWindow("Log");
				saveAs("Text", myDir+"Object centre"+mm+"mm"+mAs+"mAs");
				selectWindow("Log");
				run("Close");
			}
		}
		//CALCULATION FOR RECONSTRUCTED IMAGES
		if (ImageType_array[k]==3) //RECONSTRUCTED
		{
			nThick = nSlices;
			
			/*
			BackGround_Chest = Array.slice(Hor_Al_Profile, 0, 25);
			Array.getStatistics(BackGround_Chest, min, max, mean, std);
			BackGround_Chest_mean = mean;
			BackGround_Chest_std = std;
			
			BackGround_Nipple = Array.slice(Hor_Al_Profile, Hor_Al_Profile.length-25);
			Array.getStatistics(BackGround_Nipple, min, max, mean, std);
			BackGround_Nipple_mean = mean;
			BackGround_Nipple_std = std;
			BackGround_mean = (BackGround_Chest_mean + BackGround_Nipple_mean)/2;
			BackGround_mean_array[k] = BackGround_mean;
			BackGround_std = (BackGround_Chest_std + BackGround_Nipple_std)/2;
			BackGround_std_array[k] = BackGround_std;
			BackGround_SNR = BackGround_mean/BackGround_std;
			BackGround_SNR_array[k] = BackGround_SNR;
			
			Al = Array.slice(Hor_Al_Profile, Hor_Al_Profile.length/2-25, Hor_Al_Profile.length/2+25);
			Array.getStatistics(Al, min, max, mean, std);
			Al_mean = mean ;
			Al_mean_array[k] = Al_mean;
			Al_std = std;
			Al_std_array[k] = Al_std;
			Al_SNR = Al_mean/Al_std;
			Al_SNR_array[k] = Al_SNR;
			
			Contrast_rec = (Al_mean - BackGround_mean)/BackGround_mean;
			Contrast_rec_array[k] = Contrast_rec;
			CNR_rec = (Al_mean - BackGround_mean)/sqrt((BackGround_std*BackGround_std + Al_std*Al_std)/2);
			CNR_rec_array [k] = CNR_rec;
			*/
			

			//Homogeneity CALCULATION FOR RECONSTRUCTED IMAGES
			selectWindow(filename);
			if(LSF_cal == 1)	
			{
				setSlice(focused_slice_real);     //selection of the focused slice
			}
			else setSlice(focused_slice);
			DeltaImageBorders = floor(40/0.085);
			Homogeneity_ROI_Size = 235;
			//ROI 1 Homogeneity
			makeRectangle(DeltaImageBorders, DeltaImageBorders, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "ROI1_Homogeneity");
			run("Enhance Contrast", "saturated=0.35");
			getRawStatistics(area, mean, min, max, std);
			MeanROI1[k] = mean;
			StdROI1 = std;
			SNR_ROI1[k] = mean/std;
			
			//ROI 2 Homogeneity
			makeRectangle(nWidth-DeltaImageBorders-Homogeneity_ROI_Size, DeltaImageBorders, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "ROI2_Homogeneity");
			run("Enhance Contrast", "saturated=0.35");
			getRawStatistics(area, mean, min, max, std);
			MeanROI2[k] = mean;
			StdROI2 = std;
			SNR_ROI2[k] = mean/std;
			
			//ROI 3 Homogeneity
			makeRectangle(nWidth-DeltaImageBorders-Homogeneity_ROI_Size, nHeight-DeltaImageBorders-Homogeneity_ROI_Size, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "ROI3_Homogeneity");
			run("Enhance Contrast", "saturated=0.35");
			getRawStatistics(area, mean, min, max, std);
			MeanROI3[k] = mean;
			StdROI3 = std;
			SNR_ROI3[k] = mean/std;
			
			//ROI 4 Homogeneity
			makeRectangle(DeltaImageBorders, nHeight-DeltaImageBorders-Homogeneity_ROI_Size, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "ROI4_Homogeneity");
			run("Enhance Contrast", "saturated=0.35");
			getRawStatistics(area, mean, min, max, std);
			MeanROI4[k] = mean;
			StdROI4 = std;
			SNR_ROI4[k] = mean/std;
			
			//ROI 5 Homogeneity
			if(0.5*nWidth+0.5*Homogeneity_ROI_Size < Xobject-Object_lenght_X)
			{
				makeRectangle(0.5*nWidth-0.5*Homogeneity_ROI_Size, 0.5*nHeight-0.5*Homogeneity_ROI_Size, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI5_Homogeneity");
			}
			else
			{
				makeRectangle(Xobject-Object_lenght_X-Homogeneity_ROI_Size, 0.5*nHeight-0.5*Homogeneity_ROI_Size, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI5_Homogeneity");
			}
			run("Enhance Contrast", "saturated=0.35");
			getRawStatistics(area, mean, min, max, std);
			MeanROI5[k] = mean;
			StdROI5 = std;
			SNR_ROI5[k] = mean/std;
			
			//Percentage Deviation
			MeanTot[k] = (MeanROI1[k] + MeanROI2[k] + MeanROI3[k] + MeanROI4[k] + MeanROI5[k])/5;
			SNRTot[k] = (SNR_ROI1[k] + SNR_ROI2[k] + SNR_ROI3[k] + SNR_ROI4[k] + SNR_ROI5[k])/5;
			
			Mean_Dev1[k] = abs(1-(MeanROI1[k])/(MeanTot[k]));
			Mean_Dev2[k] = abs(1-(MeanROI2[k])/(MeanTot[k]));
			Mean_Dev3[k] = abs(1-(MeanROI3[k])/(MeanTot[k]));
			Mean_Dev4[k] = abs(1-(MeanROI4[k])/(MeanTot[k]));
			Mean_Dev5[k] = abs(1-(MeanROI5[k])/(MeanTot[k]));
		
			Mean_Dev_array[0] = Mean_Dev1[k]; 
			Mean_Dev_array[1] = Mean_Dev2[k];
			Mean_Dev_array[2] = Mean_Dev3[k]; 
			Mean_Dev_array[3] = Mean_Dev4[k]; 
			Mean_Dev_array[4] = Mean_Dev5[k];
			
			Array.getStatistics(Mean_Dev_array, min, max, mean, stdDev);
			Max_Mean_Dev[k] = max;
			
			SNR_Dev1[k] = abs(1-(SNR_ROI1[k])/(SNRTot[k]));
			SNR_Dev2[k] = abs(1-(SNR_ROI2[k])/(SNRTot[k]));
			SNR_Dev3[k] = abs(1-(SNR_ROI3[k])/(SNRTot[k]));
			SNR_Dev4[k] = abs(1-(SNR_ROI4[k])/(SNRTot[k]));
			SNR_Dev5[k] = abs(1-(SNR_ROI5[k])/(SNRTot[k]));
			
			SNR_Dev_array[0] = SNR_Dev1[k]; 
			SNR_Dev_array[1] = SNR_Dev2[k];
			SNR_Dev_array[2] = SNR_Dev3[k]; 
			SNR_Dev_array[3] = SNR_Dev4[k]; 
			SNR_Dev_array[4] = SNR_Dev5[k];
			
			Array.getStatistics(SNR_Dev_array, min, max, mean, stdDev);
			Max_SNR_Dev[k] = max;
			
			//NPS CALCULATION for RECO images
			NPS_factor = (pxdim*pxdim)/(NPS_ROI_size*NPS_ROI_size);
			//Noise Upper Panel
			DeltaX_Obj_Pos_Upper = Xobject;
			DeltaY_Obj_Pos_Upper = Yobject-4*NPS_ROI_size;

			for (i=0; i < nROI/4; i++)
			{
				for (j=0; j < nROI/4; j++)
				{
					selectWindow(filename);
					if(LSF_cal == 1)	
					{
						setSlice(focused_slice_real);     //selection of the focused slice
					}
					else setSlice(focused_slice);
					makeRectangle(DeltaX_Obj_Pos_Upper+i*0.25*NPS_ROI_size, DeltaY_Obj_Pos_Upper+j*0.25*NPS_ROI_size, NPS_ROI_size, NPS_ROI_size);
					roiManager("Add");
					roiManager("Select", roiManager("count")-1);
					roiManager("Rename", "NPS_Upper");
					getRawStatistics(area, mean, min, max, std);
					meanROI = mean;
					//Fourier Transform
					run("FFT Options...", "  raw");
					run("FFT");
					run("Line Width...", "line=5");
					makeLine(128, 128, 128, 256);
					Hor_Profile = getProfile();
					for (z=0; z<Hor_Profile.length; z++)
					{
						Hor_Profile[z] = NPS_factor*(Hor_Profile[z])/(meanROI*meanROI);
						buffer_hor[z]=buffer_hor[z]+Hor_Profile[z];
					}
					makeLine(128, 128, 256, 128);
					Ver_Profile = getProfile();
					for (z=0; z<Ver_Profile.length; z++)
					{
						Ver_Profile[z] = NPS_factor*(Ver_Profile[z])/(meanROI*meanROI);
						buffer_ver[z]=buffer_ver[z]+Ver_Profile[z];
					}
					run("Close");
				}
			}
			
			//Upper Panel
			//Calculation and printing of the average horizontal profile

			for (z=0; z<Hor_Profile.length; z++)
			{
				buffer_hor[z] = (buffer_hor[z])/nROI;
				hor_upper[z] = buffer_hor[z];
				if(dbg_val == 1) print(buffer_hor[z]);
			}
			
			if(dbg_val == 1)
			{
	
				print(d2s(mm,2));
				print(d2s(NPS_factor,-2));
				print(d2s(meanROI,2));
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Horizontal_Upper_Panel_RECO"+mm+"mm"+mAs_tot+"mAstot."+focused_slice_real+"slice");
				selectWindow("Log");
				run("Close");
			}

			//Calculation and printing of the average vertical profile

			for (z=0; z<Ver_Profile.length; z++)
			{
				buffer_ver[z] = (buffer_ver[z])/nROI;
				ver_upper[z] = buffer_ver[z];
				if(dbg_val == 1) print(buffer_ver[z]);
			}
			if(dbg_val == 1)
			{
				print(d2s(mm,2));
				print(d2s(NPS_factor,-2));
				print(d2s(meanROI,2));
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Vertical_Upper_Panel_RECO"+mm+"mm"+mAs_tot+"mAstot."+focused_slice_real+"slice");
				selectWindow("Log");
				run("Close");
			}
			//Calculation and printing of the average radial profile Upper Panel

			for (z=0; z<Ver_Profile.length; z++)
			{
				buffer_rad[z] = sqrt(buffer_ver[z]*buffer_ver[z]+buffer_hor[z]*buffer_hor[z]);
				rad_upper[z] = buffer_rad[z];
				if(dbg_val == 1) print(buffer_rad[z]);
			}
			
			if(dbg_val == 1)
			{
				print(d2s(mm,2));
				print(d2s(NPS_factor,-2));
				print(d2s(meanROI,2));
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Radial_Upper_Panel_RECO"+mm+"mm"+mAs_tot+"mAstot."+focused_slice_real+"slice");
				selectWindow("Log");
				run("Close");
			}
			//Calculation of Peak and Ground values in NPS Horizontal profile of Upper Panel;
			NPS_Ground_Upper = 0;
			for(j=40; j<50; j++)
			{
				NPS_Ground_Upper += buffer_hor[j];
			}
			NPS_Ground_Upper = NPS_Ground_Upper/10;

			//Noise Lower Panel
			DeltaX_Obj_Pos_Upper = Xobject;
			DeltaY_Obj_Pos_Upper = Yobject+2*NPS_ROI_size;

			for (i=0; i < nROI/4; i++)
			{
				for (j=0; j < nROI/4; j++)
				{
					selectWindow(filename);
					if(LSF_cal == 1)	
					{
						setSlice(focused_slice_real);     //selection of the focused slice
					}
					else setSlice(focused_slice);
					makeRectangle(DeltaX_Obj_Pos_Upper+i*0.25*NPS_ROI_size, DeltaY_Obj_Pos_Upper+j*0.25*NPS_ROI_size, NPS_ROI_size, NPS_ROI_size);
					roiManager("Add");
					roiManager("Select", roiManager("count")-1);
					roiManager("Rename", "NPS_Lower");
					getRawStatistics(area, mean, min, max, std);
					meanROI = mean;
					//Fourier Transform
					run("FFT Options...", "  raw");
					run("FFT");
					run("Line Width...", "line=5");
					makeLine(128, 128, 128, 256);
					Hor_Profile = getProfile();
					for (z=0; z<Hor_Profile.length; z++)
					{
						Hor_Profile[z] = NPS_factor*(Hor_Profile[z])/(meanROI*meanROI);
						buffer_hor[z]=buffer_hor[z]+Hor_Profile[z];
					}
					makeLine(128, 128, 256, 128);
					Ver_Profile = getProfile();
					for (z=0; z<Ver_Profile.length; z++)
					{
						Ver_Profile[z] = NPS_factor*(Ver_Profile[z])/(meanROI*meanROI);
						buffer_ver[z]=buffer_ver[z]+Ver_Profile[z];
					}
					run("Close");
				}
			}
			
			//Lower Panel
			//Calculation and printing of the average Horizontal profile Lower Panel RECO images
			for (z=0; z<Hor_Profile.length; z++)
			{
				buffer_hor[z] = (buffer_hor[z])/nROI;
				hor_lower[z] = buffer_hor[z];
				if(dbg_val == 1) print(buffer_hor[z]);
			}
			if(dbg_val == 1)
			{
				print("NPS_factor:\t"+d2s(NPS_factor,-2));
				print("Mean ROI:\t"+d2s(meanROI,2));
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Horizontal_Lower_Panel_RECO"+mm+"mm"+mAs_tot+"mAstot."+focused_slice_real+"slice");
				selectWindow("Log");
				run("Close");
			}
			//Calculation and printing of the average Vertical profile Lower Panel for RECO images
			for (z=0; z<Ver_Profile.length; z++)
			{
				buffer_ver[z] = (buffer_ver[z])/nROI;
				ver_lower[z] = buffer_ver[z];
				if(dbg_val == 1) print(buffer_ver[z]);
			}
			
			if(dbg_val == 1)
			{
				print(d2s(mm,2));
				print(d2s(NPS_factor,-2));
				print(d2s(meanROI,2));
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Vertical_Lower_Panel_RECO"+mm+"mm"+mAs_tot+"mAstot."+focused_slice_real+"slice");
				selectWindow("Log");
				run("Close");
			}
			//Calculation and printing of the average radial profile Lower Panel
			for (z=0; z<Ver_Profile.length; z++)
			{
				buffer_rad[z] = sqrt(buffer_ver[z]*buffer_ver[z]+buffer_hor[z]*buffer_hor[z]);
				rad_lower[z] = buffer_rad[z];
				if(dbg_val == 1) print(buffer_rad[z]);
			}
			if(dbg_val == 1)
			{
				print(d2s(mm,2));
				print(d2s(NPS_factor,-2));
				print(d2s(meanROI,2));
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Radial_Lower_Panel_RECO"+mm+"mm"+mAs_tot+"mAstot."+focused_slice_real+"slice");
				selectWindow("Log");
				run("Close");
			}
			//Calculation of Ground values in NPS Horizontal profile of Upper Panel for RECO images;
			NPS_Ground_Lower = 0;
			for(j=40; j<50; j++)
			{
				NPS_Ground_Lower += buffer_hor[j];
			}
			NPS_Ground_Lower = NPS_Ground_Lower/10;
			
			//Average value of the NPS Ground
			//NPS RADIAL PROFILE
			for (z=0; z<Ver_Profile.length; z++)
			{
				NPS_rad[z] = (rad_lower[z] + rad_upper[z])/2;
				print(NPS_rad[z]);
			}
			if(context == 2) //TOMO
			{
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Radial_TOMO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
				selectWindow("Log");
				run("Close");
			}
			if(context == 3) //COMBO
			{
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Radial_COMBO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
				selectWindow("Log");
				run("Close");
			}
			//NPS HORIZONTAL PROFILE
			for (z=0; z<Ver_Profile.length; z++)
			{
				NPS_hor[z] = (hor_lower[z] + hor_upper[z])/2;
				print(NPS_hor[z]);
			}
			if(context == 2) //TOMO
			{
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Horizontal_TOMO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
				selectWindow("Log");
				run("Close");
			}
			if(context == 3) //COMBO
			{
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Horizontal_COMBO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
				selectWindow("Log");
				run("Close");
			}
			//NPS VERTICAL PROFILE
			for (z=0; z<Ver_Profile.length; z++)
			{
				NPS_ver[z] = (ver_lower[z] + ver_upper[z])/2;
				print(NPS_ver[z]);
			}
			if(context == 2) //TOMO
			{
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Vertical_TOMO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
				selectWindow("Log");
				run("Close");
			}
			if(context == 3) //COMBO
			{
				selectWindow("Log");
				saveAs("Text", myDir+"NPS_Vertical_COMBO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
				selectWindow("Log");
				run("Close");
			}
			NPS_Ground_array[k] = (NPS_Ground_Upper+NPS_Ground_Lower)/2;

			//ASF Calculation
			Object_mean_asf_slice = newArray(nThick);
			Object_mean_asf_slice_sigma = newArray(nThick);
			
			ROI_1_ground_asf_slice = newArray(nThick);
			ROI_1_ground_asf_slice_sigma = newArray(nThick);
			
			ROI_2_ground_asf_slice = newArray(nThick);
			ROI_2_ground_asf_slice_sigma = newArray(nThick);
			
			Contrast_asf_slice = newArray(nThick);
			CNR_asf_slice = newArray(nThick);
			
			asf_index = newArray(nThick);
			
			DeltaPos = floor(9/0.085);
			print("Mean Al\t s Al\t ROI1\t s ROI1\t ROI2\t s ROI2\t Contrast\t CNR");
			for (s=1; s<=nThick; s++)
			{
				selectWindow(filename);
				setSlice(s);
				//Object ROI statistic
				makeRectangle(Xobject-Object_ROIsize/2, Yobject-Object_ROIsize/2, Object_ROIsize, Object_ROIsize);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_ASF_Al");
				getRawStatistics(area, mean, min, max, std);
				Object_mean_asf = mean;
				ROI_Al_ground_asf_sigma = std;
				Object_mean_asf_slice[s-1] = Object_mean_asf;
				Object_mean_asf_slice_sigma[s-1] = ROI_Al_ground_asf_sigma;
				
				makeRectangle(Xobject+DeltaPos+0.5*Object_ROIsize, Yobject-0.5*Object_ROIsize, Object_ROIsize, Object_ROIsize);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_ASF_Ground1");
				getRawStatistics(area, mean, min, max, std);
				ROI_1_ground_asf = mean;
				ROI_1_ground_asf_sigma = std;
				ROI_1_ground_asf_slice[s-1] = ROI_1_ground_asf;
				ROI_1_ground_asf_slice_sigma[s-1] = ROI_1_ground_asf_sigma;
				
				makeRectangle(Xobject-DeltaPos-2*Object_ROIsize, Yobject-0.5*Object_ROIsize, Object_ROIsize, Object_ROIsize);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_ASF_Ground2");
				getRawStatistics(area, mean, min, max, std);
				ROI_2_ground_asf = mean;
				ROI_2_ground_asf_sigma = std;
				ROI_2_ground_asf_slice[s-1] =ROI_2_ground_asf;
				ROI_2_ground_asf_slice_sigma[s-1] = ROI_2_ground_asf_sigma;
				
				Contrast_asf = (Object_mean_asf - ((ROI_1_ground_asf + ROI_2_ground_asf)/2))/((ROI_1_ground_asf + ROI_2_ground_asf)/2);
				Delta_signal_asf = (Object_mean_asf - ((ROI_1_ground_asf + ROI_2_ground_asf)/2));
				Noise_Ground_asf = (ROI_1_ground_asf_sigma+ROI_2_ground_asf_sigma)/2;
				Noise_Al_asf = ROI_Al_ground_asf_sigma;
				CNR_asf = Delta_signal_asf/sqrt((Noise_Ground_asf*Noise_Ground_asf + Noise_Al_asf*Noise_Al_asf)/2);
				Contrast_asf_slice[s-1] = Contrast_asf;
				CNR_asf_slice[s-1] = CNR_asf;
				asf_index[s-1] = s;
				print(Object_mean_asf_slice[s-1]+"\t"+Object_mean_asf_slice_sigma[s-1]+"\t"+ROI_1_ground_asf_slice[s-1]+"\t"+ROI_1_ground_asf_slice_sigma[s-1]+"\t"+ROI_2_ground_asf_slice[s-1]+"\t"+ROI_2_ground_asf_slice_sigma[s-1]+"\t"+d2s(Contrast_asf_slice[s-1],3)+"\t"+d2s(CNR_asf_slice[s-1],3));
			}
			
			if(LSF_cal == 1)
			{
				BackGround_mean_array[k] = (ROI_1_ground_asf_slice[focused_slice_real-1] + ROI_2_ground_asf_slice[focused_slice_real-1])/2;
				BackGround_std_array[k] = (ROI_1_ground_asf_slice_sigma[focused_slice_real-1] + ROI_2_ground_asf_slice_sigma[focused_slice_real-1])/2;
				BackGround_SNR_array[k] = (BackGround_mean_array[k])/(BackGround_std_array[k]);
				
				Al_mean_array[k] = Object_mean_asf_slice[focused_slice_real-1];
				Al_std_array[k] = Object_mean_asf_slice_sigma[focused_slice_real-1];
				Al_SNR_array[k] = (Al_mean_array[k])/(Al_std_array[k]);
				
				Contrast_rec_array[k] = Contrast_asf_slice[focused_slice_real-1];
				CNR_asf_focused_slice[k] = CNR_asf_slice[focused_slice_real-1];
			}
			else
			{
				BackGround_mean_array[k] = (ROI_1_ground_asf_slice[focused_slice-1] + ROI_2_ground_asf_slice[focused_slice-1])/2;
				BackGround_std_array[k] = (ROI_1_ground_asf_slice_sigma[focused_slice-1] + ROI_2_ground_asf_slice_sigma[focused_slice-1])/2;
				BackGround_SNR_array[k] = (BackGround_mean_array[k])/(BackGround_std_array[k]);
				
				Al_mean_array[k] = Object_mean_asf_slice[focused_slice-1];
				Al_std_array[k] = Object_mean_asf_slice_sigma[focused_slice-1];
				Al_SNR_array[k] = (Al_mean_array[k])/(Al_std_array[k]);
				
				Contrast_rec_array[k] = Contrast_asf_slice[focused_slice-1];
				CNR_asf_focused_slice[k] = CNR_asf_slice[focused_slice-1];
			}
			Fit.doFit("gaussian", asf_index,Contrast_asf_slice);
			a = Fit.p(0);
			b = Fit.p(1);
			mu = Fit.p(2);
			sigma = Fit.p(3);
			sigma_asf[k] = sigma;
			if(context == 2) //TOMO
			{
				print(a);
				print(b);
				print(mu);
				print(sigma);
				selectWindow("Log");
				saveAs("Text", myDir+"ASF_TOMO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
				selectWindow("Log");
				run("Close");
			}
			if(context == 3) //COMBO
			{
				print(a);
				print(b);
				print(mu);
				print(sigma);
				selectWindow("Log");
				saveAs("Text", myDir+"ASF_COMBO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
				selectWindow("Log");
				run("Close");
			//Test result
			}
			if (Filter==" SILVER" && BackGround_SNR_array[k] > 50)
			{
				result_SNR = "Pass";
			}
			if (Filter==" RHODIUM " && BackGround_SNR_array[k] > 50)
			{
				result_SNR = "Pass";
			}
			result_SNR_array[k] = result_SNR;
			
			//LSF Calculation for RECO images
			if(LSF_cal ==1)
			{
				//Detection of the tungsten wire
				LSF_Offset_X = 35/pxdim; //Position of the tungsten wire with respect to the Al object in mm on X direction
				LSF_Offset_Y = 10/pxdim; //Position of the tungsten wire with respect to the Al object in mm on Y direction
				Wire_detection_ROI = 200; //ROI dimension for tungsten wire detection 
				
				selectWindow(filename);
				setSlice(focused_slice_real);
				makeRectangle(Xobject+LSF_Offset_X-0.5*Wire_detection_ROI, Yobject-LSF_Offset_Y-0.5*Wire_detection_ROI, Wire_detection_ROI, Wire_detection_ROI);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_LSF");
				run("Copy");
				newImage("copy_image_for_wire_detection", "16-bit Black", Wire_detection_ROI, Wire_detection_ROI, 1);
				run("Paste");
				run("Enhance Contrast", "saturated=0.35");
				run("Select All");
				getRawStatistics(area, mean, min, max, std);
				Mean_Wire_detection_ROI = mean;
				run("Subtract...", "value=Mean_Wire_detection_ROI");
				run("Enhance Contrast", "saturated=0.35");

				//Detection of the peak and coordinates calculation
				x = newArray(Wire_detection_ROI);
				Peak = newArray(Wire_detection_ROI);
				y = newArray(Wire_detection_ROI); 
				for(i = 0; i < Wire_detection_ROI; i++)
				{
					x[i] = Wire_detection_ROI-i;
					makeLine(i,0,i,Wire_detection_ROI);
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
				if (dbg_val ==1 )
				{
					selectWindow("Log");
					saveAs("Text", myDir+"Coordinates wire Left-Right direction"+mm+"mm");
					selectWindow("Log");
					run("Close");
				}
				Fit.doFit("straight line", x, y);
				m = Fit.p(1);
				q = Fit.p(0);
				R_square_lin = Fit.rSquared;
				angle_wire = -atan(Fit.p(1));
				angle_wire_degree = atan(Fit.p(1))*180/PI;
				selectWindow("copy_image_for_wire_detection");
				run("Close");
				
				//Perpendicular Profiles of the wire 
				selectWindow(filename);
				setSlice(focused_slice_real);
				run("Enhance Contrast", "saturated=0.35");
				makeRectangle(Xobject+LSF_Offset_X-0.5*Wire_detection_ROI, Yobject-LSF_Offset_Y-0.5*Wire_detection_ROI, Wire_detection_ROI, Wire_detection_ROI);
				run("Copy");
				newImage("copy_image_for_LSF_calculation", "16-bit Black", Wire_detection_ROI, Wire_detection_ROI, 1);
				run("Paste");
				run("Enhance Contrast", "saturated=0.35");
				selectWindow("copy_image_for_LSF_calculation");
				
				LSF_length = 128; //length of the wired profile;
				a_array = newArray(Wire_detection_ROI-20); //array containing the value of the gaussian fit parameters a
				b_array = newArray(Wire_detection_ROI-20); //array containing the value of the gaussian fit parameters b
				mu_array = newArray(Wire_detection_ROI-20); //array containing the value of the gaussian fit parameters mu
				sigma_array = newArray(Wire_detection_ROI-20); //array containing the value of the gaussian fit parameters sigma
				R_square_array = newArray(Wire_detection_ROI-20); //array containing the value of the gaussian fit parameters R_square
				lsf_left_right_sum = newArray(LSF_length);
				lsf_left_right_media = newArray(LSF_length);
				X_lsf = newArray(LSF_length);
				for(i = 10; i < Wire_detection_ROI-10; i++)
				{
					
					x[i] = i;
					XStart = x[i]+LSF_length/2*sin(angle_wire);
					YStart = y[i]-LSF_length/2*cos(angle_wire);
					XEnd = x[i]-LSF_length/2*sin(angle_wire);
					YEnd = y[i]+LSF_length/2*cos(angle_wire);
					run("Line Width...", "line=2");
					makeLine(XStart,YStart,XEnd,YEnd);
					roiManager("Add");
					roiManager("Select", roiManager("count")-1);
					roiManager("Rename", "ROI_LSF_profile"+i-9);
					lsf_left_right = getProfile();
					for (u=0; u < lsf_left_right_sum.length; u++)
					{
						lsf_left_right_sum[u] += lsf_left_right[u];
					}
				
					for(j=0; j < lsf_left_right.length; j++)
					{
					
						X_lsf[j] = (j+1)*pxdim; //mm
					}
					
					Fit.doFit("gaussian", X_lsf,lsf_left_right);
					R_square = Fit.rSquared;
					a = Fit.p(0);
					b = Fit.p(1);
					mu = Fit.p(2);
					sigma = Fit.p(3);
					a_array[i-10] = a;
					b_array[i-10] = b;
					mu_array[i-10] = mu;
					sigma_array[i-10] = abs(sigma);
					R_square_array[i-10] = R_square;
					
					sigma_new = newArray(Wire_detection_ROI-20);
					a_new = newArray(Wire_detection_ROI-20);
					b_new = newArray(Wire_detection_ROI-20);
					mu_new = newArray(Wire_detection_ROI-20);
					new_count = 0;
					//Exclude the values for R^2 <0.4
					for(e=0; e < sigma_new.length; e++)
					{
						if(R_square_array[e] > 0.4)
						{
							a_new[new_count] = a_array[e];
							b_new[new_count] = b_array[e];
							mu_new[new_count] = mu_array[e];
							sigma_new[new_count] = sigma_array[e];
							new_count++;
						}
					}
				}
				
				//LSF average profile calculation and printing
				for (u=0; u < lsf_left_right.length; u++)
				{
					lsf_left_right_media[u] = lsf_left_right_sum[u]/(Wire_detection_ROI-20);
					if(dbg_val ==1) print(d2s(X_lsf[u],4)+"\t"+d2s(lsf_left_right_media[u],2));
				}
				Fit.doFit("gaussian", X_lsf,lsf_left_right_media);
				R_square_gauss_mean = Fit.rSquared;
				a_gauss_mean = Fit.p(0);
				b_gauss_mean = Fit.p(1);
				mu_gauss_mean = Fit.p(2);
				sigma_gauss_mean = Fit.p(3);
				
				if(dbg_val ==1)
				{
					print(d2s(mm,2));
					print(d2s(mAs,2));
					print(d2s(a_gauss_mean,2));
					print(d2s(b_gauss_mean,2));
					print(d2s(mu_gauss_mean,2));
					print(d2s(sigma_gauss_mean,3));
					print(d2s(R_square_gauss_mean,3));
					if(context == 2 && dbg_val ==1) //TOMO
					{
						selectWindow("Log");
						saveAs("Text", myDir+"LSF_mean_TOMO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
						selectWindow("Log");
						run("Close");
					}
					if(context == 3 && dbg_val ==1) //COMBO
					{
						selectWindow("Log");
						saveAs("Text", myDir+"LSF_mean_COMBO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
						selectWindow("Log");
						run("Close");
					}
				}
				a_stat = Array.trim(a_new, new_count);
				b_stat = Array.trim(b_new, new_count);
				mu_stat = Array.trim(mu_new, new_count);
				sigma_stat = Array.trim(sigma_new, new_count);
				
				Array.getStatistics(a_stat, min, max, mean, std);
				a_mean = mean;
				a_std = std;
				
				Array.getStatistics(b_stat, min, max, mean, std);
				b_mean = mean;
				b_std = std;
				
				Array.getStatistics(mu_stat, min, max, mean, std);
				mu_mean = mean;
				mu_std = std;

				Array.getStatistics(sigma_stat, min, max, mean, std);
				sigma_mean = mean;
				sigma_std = std;
				
				stat_data_number = sigma_stat.length;
				if (dbg_val == 1)
				{
					print("Linear Fit parameters:\t");
					print("m:\t"+d2s(m,3));
					print("q:\t"+d2s(q,3));
					print("R2:\t"+d2s(R_square_lin,3));
					print("angle (degree):\t"+d2s(angle_wire_degree,3));
					
					print("Gaussian Fit parameter:");
					print("Used"+stat_data_number+"profiles");
					print("a:\t"+d2s(a_mean,0)+"\t"+d2s(a_std,0));
					print("b:\t"+d2s(b_mean,-0)+"\t"+d2s(b_std,0));
					print("mu:\t"+d2s(mu_mean,3)+"\t"+d2s(mu_std,3));
					print("sigma mean:\t"+d2s(sigma_mean,3)+"\t"+d2s(sigma_std,3));
					
					print("Fit parameter for Average Gaussian curve:");
					print("a:\t"+d2s(a_gauss_mean,0));
					print("b:\t"+d2s(b_gauss_mean,-0));
					print("mu:\t"+d2s(mu_gauss_mean,3));
					print("sigma mean:\t"+d2s(sigma_gauss_mean,3));
					print("R2:\t"+d2s(R_square_gauss_mean,3));
				}
				scale = 2*sqrt(2*log(2));
				lsf_array[k] = abs(sigma_gauss_mean*scale);
				
				if (context == 2 && dbg_val == 1)//TOMO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"LSF_Result_Left-Right_TOMO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
					selectWindow("Log");
					run("Close");
				}
				if (context == 3 && dbg_val == 1)//COMBO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"LSF_Result_Left-Right_COMBO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
					selectWindow("Log");
					run("Close");
				}
				
				//MTF Calculation for TOMO or COMBO RECO images
				//Pre-sampling LSF variables
				
				startpt = newArray(50, 75, 100, 125); //Starting point for the 4 LSF profiles
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
				xoff = newArray(1/sin(abs(angle_wire)));
				xOC = newArray(lsf.length*xoff.length);
				lsfOC = newArray(lsf.length*xoff.length);
				lsf_chest_nipple_sum = newArray(129);
				lsf_chest_nipple_media = newArray(129);
				wire_angle_array[k] = angle_wire_degree;
				//Calculation of the Pre-sampling LSF
				
				selectWindow(filename);
				setSlice(focused_slice_real);
				run("Enhance Contrast", "saturated=0.35");
				makeRectangle(Xobject+LSF_Offset_X-0.5*Wire_detection_ROI, Yobject-LSF_Offset_Y-0.5*Wire_detection_ROI, Wire_detection_ROI, Wire_detection_ROI);
				run("Copy");
				newImage("copy_image_for_Pre_sampling_LSF_calculation", "16-bit Black", Wire_detection_ROI, Wire_detection_ROI, 1);
				run("Paste");
				run("Enhance Contrast", "saturated=0.35");
				selectWindow("copy_image_for_Pre_sampling_LSF_calculation");
				
				
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

					if(context == 2 && dbg_val == 1)
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
						saveAs("Text", myDir+"Pre_sampling_LSF_raw_TOMO_REC"+z+1+"lsf"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
						selectWindow("Log");
						run("Close");
					}
					if(context == 3 && dbg_val == 1)
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
						saveAs("Text", myDir+"Pre_sampling_LSF_raw_COMBO_REC"+z+1+"lsf"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
						selectWindow("Log");
						run("Close");
					}
					LSF_OC_stat_a[z] = a;
					LSF_OC_stat_b[z] = b;
					LSF_OC_stat_c[z] = c;
					LSF_OC_stat_d[z] = abs(d);
					LSF_OC_stat_FHWM[z] = abs(scale*d);					//vector containing the sigma value of the pre-sampling lsf
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

					//0 value for the pixel under 20 % of the peak and gaussian fit value for the tails
					for(j=0; j<lsfOC.length; j++)
					{
						if(lsfOC[j] <= 0.1)
						{
							lsfOC[j] = 0; //Put at 0 the ground value
						}
						if(lsfOC[j] > 0.1 && lsfOC[j] < 0.25)
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
					if(xOC_sort.length >= 2048)
					{
						slice_number = 1024;
					}
					if(xOC_sort.length >= 1024 && xOC_sort.length <= 2048)
					{
						slice_number = 512;
					}
					if(xOC_sort.length >= 512 && xOC_sort.length <= 1024)
					{
						slice_number = 256;
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
					if(context == 2 && dbg_val == 1)
					{
						print("a\t"+a);
						print("mu\t"+b);
						print("sigma\t"+c);
						selectWindow("Log");
						saveAs("Text", myDir+"Pre_sampling_LSF_normalized_buffer_TOMO"+z+1+"lsf"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
						selectWindow("Log");
						run("Close");
					}
					if(context == 3 && dbg_val == 1)
					{
						print("a\t"+a);
						print("mu\t"+b);
						print("sigma\t"+c);
						selectWindow("Log");
						saveAs("Text", myDir+"Pre_sampling_LSF_normalized_buffer_COMBO"+z+1+"lsf"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
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
					if(context == 2 && dbg_val == 1)
					{
						for (i =0; i< modulo_norm.length; i++)
						{
							print(modulo_norm[i]);
						}
						selectWindow("Log");
						saveAs("Text", myDir+"MTF_TOMO_RECO"+(z+1)+"_"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
						selectWindow("Log");
						run("Close");
					}
					if(context == 3 && dbg_val == 1)
					{
						for (i =0; i< modulo_norm.length; i++)
						{
							print(modulo_norm[i]);
						}
						selectWindow("Log");
						saveAs("Text", myDir+"MTF_COMBO_RECO"+(z+1)+"_"+mm+"mm"+d2s(mAs,0)+"mAs");
						selectWindow("Log");
						run("Close");
					}
				}
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
				if(context == 2)
				{
					for(i=0; i<lsf_buffer.length; i++)
					{
						print(d2s(x_to_print[i],3)+"\t"+d2s(lsf_to_print[i],2)+"\t"+d2s(a_lsf_to_print*exp(-((x_to_print[i]-b_lsf_to_print)*(x_to_print[i]-b_lsf_to_print))/(2*c_lsf_to_print*c_lsf_to_print)),2));
					}
					print("a\t"+d2s(a_lsf_to_print,2));
					print("b\t"+d2s(b_lsf_to_print,2));
					print("c\t"+d2s(c_lsf_to_print,2));
					selectWindow("Log");
					saveAs("Text", myDir+"LSF_TOMO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3)
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
					saveAs("Text", myDir+"LSF_COMBO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
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
				if(context == 2)
				{
					selectWindow("Log");
					saveAs("Text", myDir+"MTF_TOMO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3)
				{
					selectWindow("Log");
					saveAs("Text", myDir+"MTF_COMBO_RECO"+mm+"mm"+d2s(mAs_tot,0)+"mAstot.");
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
				if(context == 2)
				{
					selectWindow("Log");
					saveAs("Text", myDir+"Pre_sampling_LSF_Fit_Results_TOMO_RECO"+mm+"mm"+mAs_tot+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3)
				{
					selectWindow("Log");
					saveAs("Text", myDir+"Pre_sampling_LSF_Fit_Results_COMBO_RECO"+mm+"mm"+mAs_tot+"mAs");
					selectWindow("Log");
					run("Close");
				}
			
			}	
			
			//Test Results
			
			//Homogeneity Result and ROISet saving
			if ((mm >35) && (mm <45))
			{
				if (Max_Mean_Dev[k] < 0.1)
				{
					result_Homogeneity_Mean = "Pass";
				}
				if (Max_Mean_Dev[k] < 0.35)
				{
					result_Homogeneity_SNR = "Pass";
				}
			}
			if(context == 2)//TOMO
			{
				if(dbg_val ==1) roiManager("Save", myDir+"RoiSet_TOMO_RECO"+mm+"mm"+mAs_tot+"mAs"+angle+"degree"+".zip");
				//roiManager("Delete");
			}
			if(context == 3)//COMBO
			{
				if(dbg_val ==1) roiManager("Save", myDir+"RoiSet_COMBO_RECO"+mm+"mm"+mAs_tot+"mAs"+angle+"degree"+".zip");
				//roiManager("Delete");
			}

		}
		else //for mammo and projection images
		{
			//Object ROI statistic
			selectWindow(filename);
			makeRectangle(Xobject-Object_ROIsize/2, Yobject-Object_ROIsize/2, Object_ROIsize, Object_ROIsize);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "ROI_Al");
			run("Enhance Contrast", "saturated=0.35");
			getRawStatistics(area, mean, min, max, std, histogram);
			StdDetail = std;
			StdDetail_array[k] = StdDetail;
			VarDetail = std*std;
			MeanDetail = mean;
			MeanDetail_array[k] = MeanDetail;
			
			//Background statistic
			DeltaPos = floor(9/0.085);
			//ROI1 (Right with respect to the Al object)
			makeRectangle(Xobject+DeltaPos+0.5*Object_ROIsize, Yobject-0.5*Object_ROIsize, Object_ROIsize, Object_ROIsize);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "ROI_CNR_1");
			run("Enhance Contrast", "saturated=0.35");
			getRawStatistics(area, mean, min, max, std, histogram);
			MeanBKG1 = mean;
			StdBKG1 = std;
			VarBKG1 = std*std;
			SNR1 = mean/std;
			
			//ROI3 (Left with respect to the Al object)
			makeRectangle(Xobject-DeltaPos-2*Object_ROIsize, Yobject-0.5*Object_ROIsize, Object_ROIsize, Object_ROIsize);
			roiManager("Add");
			roiManager("Select", roiManager("count")-1);
			roiManager("Rename", "ROI_CNR_3");
			run("Enhance Contrast", "saturated=0.35");
			getRawStatistics(area, mean, min, max, std, histogram);
			MeanBKG3 = mean;
			StdBKG3 = std;
			VarBKG3 = std*std;
			SNR3 = mean/std;
			if(	ImageType_array[k] == 1) //MAMMO
			{
				//ROI2 (Over with respect to the Al object)
				makeRectangle(Xobject-0.5*Object_ROIsize, Yobject-DeltaPos-2*Object_ROIsize, Object_ROIsize, Object_ROIsize);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_CNR_2");
				run("Enhance Contrast", "saturated=0.35");
				getRawStatistics(area, mean, min, max, std, histogram);
				MeanBKG2 = mean;
				StdBKG2 = std;
				VarBKG2 = std*std;
				SNR2 = mean/std;
				
				//ROI4 (Under with respect to the Al object)
				makeRectangle(Xobject-0.5*Object_ROIsize, Yobject+DeltaPos+Object_ROIsize, Object_ROIsize, Object_ROIsize);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_CNR_4");
				run("Enhance Contrast", "saturated=0.35");
				getRawStatistics(area, mean, min, max, std, histogram);
				MeanBKG4 = mean;
				StdBKG4 = std;
				VarBKG4 = std*std;
				SNR4 = mean/std;
			
				MeanBKG =(MeanBKG1+MeanBKG2+MeanBKG3+MeanBKG4)/4;
				MeanBKG_array[k] = MeanBKG;
				StdBKG = (StdBKG1+StdBKG2+StdBKG3+StdBKG4)/4;
				StdBKG_array[k] = StdBKG;
				VarBKG = (VarBKG1+VarBKG2+VarBKG3+VarBKG4)/4;
			}
			if(	ImageType_array[k] == 2) //PROJECTION
			{
				MeanBKG =(MeanBKG1+MeanBKG3)/2;
				MeanBKG_array[k] = MeanBKG;
				StdBKG = (StdBKG1+StdBKG3)/2;
				StdBKG_array[k] = StdBKG;
				VarBKG = (VarBKG1+VarBKG3)/2;
			}
			Contrast = abs((MeanBKG-MeanDetail)/MeanBKG);
			Contrast_array[k] = Contrast;
			Noise = sqrt((VarDetail+VarBKG)/2);
			SNR_Object=MeanDetail/StdDetail;
			SNR_Object_array[k] = SNR_Object;
			SNR_BKG=MeanBKG/StdDetail;
			SNR_BKG_array[k] = SNR_BKG;
			CNR = (MeanBKG-MeanDetail)/Noise;
			CNR_array[k] = CNR;
			
			//Homogeneity CALCULATION for only MAMMO
			if(	ImageType_array[k] == 1) //MAMMO
			{
				DeltaImageBorders = floor(30/0.085);
				Homogeneity_ROI_Size = 235;
				//ROI 1 Homogeneity
				makeRectangle(DeltaImageBorders, DeltaImageBorders, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_Homogenety_1");
				run("Enhance Contrast", "saturated=0.35");
				getRawStatistics(area, mean, min, max, std);
				MeanROI1[k] = mean;
				StdROI1 = std;
				SNR_ROI1[k] = mean/std;
				
				//ROI 2 Homogeneity
				makeRectangle(nWidth-DeltaImageBorders-Homogeneity_ROI_Size, DeltaImageBorders, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_Homogenety_2");
				run("Enhance Contrast", "saturated=0.35");
				getRawStatistics(area, mean, min, max, std);
				MeanROI2[k] = mean;
				StdROI2 = std;
				SNR_ROI2[k] = mean/std;
				
				//ROI 3 Homogeneity
				makeRectangle(nWidth-DeltaImageBorders-Homogeneity_ROI_Size, nHeight-DeltaImageBorders-Homogeneity_ROI_Size, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_Homogenety_3");
				run("Enhance Contrast", "saturated=0.35");
				getRawStatistics(area, mean, min, max, std);
				MeanROI3[k] = mean;
				StdROI3 = std;
				SNR_ROI3[k] = mean/std;
				
				//ROI 4 Homogeneity
				makeRectangle(DeltaImageBorders, nHeight-DeltaImageBorders-Homogeneity_ROI_Size, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_Homogenety_4");
				run("Enhance Contrast", "saturated=0.35");
				getRawStatistics(area, mean, min, max, std);
				MeanROI4[k] = mean;
				StdROI4 = std;
				SNR_ROI4[k] = mean/std;
				
				//ROI 5 Homogeneity
				if(0.5*nWidth+0.5*Homogeneity_ROI_Size < Xobject-Object_lenght_X)
				{
					makeRectangle(0.5*nWidth-0.5*Homogeneity_ROI_Size, 0.5*nHeight-0.5*Homogeneity_ROI_Size, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
					roiManager("Add");
					roiManager("Select", roiManager("count")-1);
					roiManager("Rename", "ROI_Homogenety_5");
				}
				else
				{
					makeRectangle(Xobject-Object_lenght_X-Homogeneity_ROI_Size, 0.5*nHeight-0.5*Homogeneity_ROI_Size, Homogeneity_ROI_Size, Homogeneity_ROI_Size);
					roiManager("Add");
					roiManager("Select", roiManager("count")-1);
					roiManager("Rename", "ROI_Homogenety_5");
				}
				run("Enhance Contrast", "saturated=0.35");
				getRawStatistics(area, mean, min, max, std);
				MeanROI5[k] = mean;
				StdROI5 = std;
				SNR_ROI5[k] = mean/std;
				
				//Percentage Deviation
				MeanTot[k] = (MeanROI1[k] + MeanROI2[k] + MeanROI3[k] + MeanROI4[k] + MeanROI5[k])/5;
				SNRTot[k] = (SNR_ROI1[k] + SNR_ROI2[k] + SNR_ROI3[k] + SNR_ROI4[k] + SNR_ROI5[k])/5;
				
				Mean_Dev1[k] = abs(1-(MeanROI1[k])/(MeanTot[k]));
				Mean_Dev2[k] = abs(1-(MeanROI2[k])/(MeanTot[k]));
				Mean_Dev3[k] = abs(1-(MeanROI3[k])/(MeanTot[k]));
				Mean_Dev4[k] = abs(1-(MeanROI4[k])/(MeanTot[k]));
				Mean_Dev5[k] = abs(1-(MeanROI5[k])/(MeanTot[k]));
				
				Mean_Dev_array[0] = Mean_Dev1[k]; 
				Mean_Dev_array[1] = Mean_Dev2[k];
				Mean_Dev_array[2] = Mean_Dev3[k]; 
				Mean_Dev_array[3] = Mean_Dev4[k]; 
				Mean_Dev_array[4] = Mean_Dev5[k];
				
				Array.getStatistics(Mean_Dev_array, min, max, mean, stdDev);
				Max_Mean_Dev[k] = max;
				
				SNR_Dev1[k] = abs(1-(SNR_ROI1[k])/(SNRTot[k]));
				SNR_Dev2[k] = abs(1-(SNR_ROI2[k])/(SNRTot[k]));
				SNR_Dev3[k] = abs(1-(SNR_ROI3[k])/(SNRTot[k]));
				SNR_Dev4[k] = abs(1-(SNR_ROI4[k])/(SNRTot[k]));
				SNR_Dev5[k] = abs(1-(SNR_ROI5[k])/(SNRTot[k]));
				
				SNR_Dev_array[0] = SNR_Dev1[k]; 
				SNR_Dev_array[1] = SNR_Dev2[k];
				SNR_Dev_array[2] = SNR_Dev3[k]; 
				SNR_Dev_array[3] = SNR_Dev4[k]; 
				SNR_Dev_array[4] = SNR_Dev5[k];
				
				Array.getStatistics(SNR_Dev_array, min, max, mean, stdDev);
				Max_SNR_Dev[k] = max;
			}
			if(ImageType_array[k] ==1) //Only for MAMMO
			//NPS CALCULATION for MAMMO images
			{
				NPS_factor = (pxdim*pxdim)/(NPS_ROI_size*NPS_ROI_size);
				
				//Noise Upper Panel
				DeltaX_Obj_Pos_Upper = Xobject;
				DeltaY_Obj_Pos_Upper = Yobject-4*NPS_ROI_size;

				for (i=0; i < nROI/4; i++)
				{
					for (j=0; j < nROI/4; j++)
					{
						selectWindow(filename);
						makeRectangle(DeltaX_Obj_Pos_Upper+i*0.25*NPS_ROI_size, DeltaY_Obj_Pos_Upper+j*0.25*NPS_ROI_size, NPS_ROI_size, NPS_ROI_size);
						roiManager("Add");
						roiManager("Select", roiManager("count")-1);
						roiManager("Rename", "ROI_NPS_Upper"+i+j);
						getRawStatistics(area, mean, min, max, std);
						meanROI = mean;
						//Fourier Transform
						run("FFT Options...", "  raw");
						run("FFT");
						run("Line Width...", "line=5");
						makeLine(128, 128, 128, 256);
						Hor_Profile = getProfile();
						for (z=0; z<Hor_Profile.length; z++)
						{
							Hor_Profile[z] = NPS_factor*(Hor_Profile[z])/(meanROI*meanROI);
							buffer_hor[z]=buffer_hor[z]+Hor_Profile[z];
						}
						makeLine(128, 128, 256, 128);
						Ver_Profile = getProfile();
						for (z=0; z<Ver_Profile.length; z++)
						{
							Ver_Profile[z] = NPS_factor*(Ver_Profile[z])/(meanROI*meanROI);
							buffer_ver[z]=buffer_ver[z]+Ver_Profile[z];
						}
						run("Close");
					}
				}
				
				//Upper Panel
				//Calculation and printing of the average horizontal profile and grid peak research

				for (z=0; z<Hor_Profile.length; z++)
				{
					buffer_hor[z] = (buffer_hor[z])/nROI;
					hor_upper[z] = buffer_hor[z];
					if(dbg_val == 1) print(buffer_hor[z]);
				}
				if(dbg_val == 1)
				{
					print(d2s(mm,2));
					print(d2s(angle,2));
					print(d2s(mAs,2));
					print(d2s(NPS_factor,-2));
					print(d2s(meanROI,2));
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Horizontal_Upper_Panel"+mm+"mm"+mAs+"mAs"+angle+"degree");
					selectWindow("Log");
					run("Close");
				}

				//Calculation and printing of the average vertical profile

				for (z=0; z<Ver_Profile.length; z++)
				{
					buffer_ver[z] = (buffer_ver[z])/nROI;
					ver_upper[z] = buffer_ver[z];
					if(dbg_val == 1) print(buffer_ver[z]);
				}
				if(dbg_val == 1)
				{
					print(d2s(mm,2));
					print(d2s(angle,2));
					print(d2s(mAs,2));
					print(d2s(NPS_factor,-2));
					print(d2s(meanROI,2));
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Vertical_Upper_Panel"+mm+"mm"+mAs+"mAs"+angle+"degree");
					selectWindow("Log");
					run("Close");
				}
				//Calculation and printing of the average radial profile Upper Panel

				for (z=0; z<Ver_Profile.length; z++)
				{
					buffer_rad[z] = sqrt(buffer_ver[z]*buffer_ver[z]+buffer_hor[z]*buffer_hor[z]);
					rad_upper[z] = buffer_rad[z];
					if(dbg_val == 1) print(buffer_rad[z]);
				}
				if(dbg_val == 1)
				{
					print(d2s(mm,2));
					print(d2s(NPS_factor,-2));
					print(d2s(meanROI,2));
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Radial_Upper_Panel_PROJECTION"+mm+"mm"+mAs+"mAs"+angle+"degree");
					selectWindow("Log");
					run("Close");
				}
				//Calculation of Peak and Ground values in NPS Horizontal profile of Upper Panel;
				Grid_Peak_Upper = buffer_hor[66];
				Grid_Peak_array[k] = Grid_Peak_Upper;
				NPS_Ground_Upper = 0;
				for(j=40; j<50; j++)
				{
					NPS_Ground_Upper += buffer_hor[j];
				}
				NPS_Ground_Upper = NPS_Ground_Upper/10;
				Peak_Ground_Ratio_Upper = Grid_Peak_Upper/NPS_Ground_Upper;

				//Noise Lower Panel
				DeltaX_Obj_Pos_Upper = Xobject;
				DeltaY_Obj_Pos_Upper = Yobject+2*NPS_ROI_size;

				for (i=0; i < nROI/4; i++)
				{
					for (j=0; j < nROI/4; j++)
					{
						selectWindow(filename);
						makeRectangle(DeltaX_Obj_Pos_Upper+i*0.25*NPS_ROI_size, DeltaY_Obj_Pos_Upper+j*0.25*NPS_ROI_size, NPS_ROI_size, NPS_ROI_size);
						roiManager("Add");
						roiManager("Select", roiManager("count")-1);
						roiManager("Rename", "ROI_NPS_Lower"+i+j);
						getRawStatistics(area, mean, min, max, std);
						meanROI = mean;
						//Fourier Transform
						run("FFT Options...", "  raw");
						run("FFT");
						run("Line Width...", "line=5");
						makeLine(128, 128, 128, 256);
						Hor_Profile = getProfile();
						for (z=0; z<Hor_Profile.length; z++)
						{
							Hor_Profile[z] = NPS_factor*(Hor_Profile[z])/(meanROI*meanROI);
							buffer_hor[z]=buffer_hor[z]+Hor_Profile[z];
						}
						makeLine(128, 128, 256, 128);
						Ver_Profile = getProfile();
						for (z=0; z<Ver_Profile.length; z++)
						{
							Ver_Profile[z] = NPS_factor*(Ver_Profile[z])/(meanROI*meanROI);
							buffer_ver[z]=buffer_ver[z]+Ver_Profile[z];
						}
						run("Close");
					}
				}
				
				//Lower Panel
				//Calculation and printing of the average horizontal profile and grid peak research

				for (z=0; z<Hor_Profile.length; z++)
				{
					buffer_hor[z] = (buffer_hor[z])/nROI;
					hor_lower[z] = buffer_hor[z];
					if(dbg_val == 1) print(buffer_hor[z]);
				}
				if(dbg_val == 1)
				{
					print(d2s(mm,2));
					print(d2s(angle,2));
					print(d2s(mAs,2));
					print(d2s(NPS_factor,-2));
					print(d2s(meanROI,2));
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Horizontal_Lower_Panel"+mm+"mm"+mAs+"mAs"+angle+"degree");
					selectWindow("Log");
					run("Close");
				}

				//Calculation and printing of the average vertical profile for Lower panel
				
				for (z=0; z<Ver_Profile.length; z++)
				{
					buffer_ver[z] = (buffer_ver[z])/nROI;
					ver_lower[z] = buffer_ver[z];
					if(dbg_val == 1) print(buffer_ver[z]);
				}
				if(dbg_val == 1)
				{
					print(d2s(mm,2));
					print(d2s(angle,2));
					print(d2s(mAs,2));
					print(d2s(NPS_factor,-2));
					print(d2s(meanROI,2));
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Vertical_Lower_Panel"+mm+"mm"+mAs+"mAs"+angle+"degree");
					selectWindow("Log");
					run("Close");
				}
				//Calculation and printing of the average radial profile Lower Panel

				for (z=0; z<Ver_Profile.length; z++)
				{
					buffer_rad[z] = sqrt(buffer_ver[z]*buffer_ver[z]+buffer_hor[z]*buffer_hor[z]);
					rad_lower[z] = buffer_rad[z];
					if(dbg_val == 1) print(buffer_rad[z]);
				}
				if(dbg_val == 1)
				{
					print(d2s(mm,2));
					print(d2s(NPS_factor,-2));
					print(d2s(meanROI,2));
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Radial_Lower_Panel_PROJECTION"+mm+"mm"+mAs+"mAs"+angle+"degree");
					selectWindow("Log");
					run("Close");
				}
				//Calculation of Peak and Ground values in NPS Horizontal profile of Lower Panel;
				Grid_Peak_Lower = buffer_hor[66];
				Grid_Peak_array[k] = Grid_Peak_Lower;
				NPS_Ground_Lower = 0;
				for(j=40; j<50; j++)
				{
					NPS_Ground_Lower += buffer_hor[j];
				}
				NPS_Ground_Lower = NPS_Ground_Lower/10;
				Peak_Ground_Ratio_Lower = Grid_Peak_Lower/NPS_Ground_Lower;
				NPS_Ground_array[k] = (NPS_Ground_Upper+NPS_Ground_Lower)/2;
				grid_visibility_array[k] = (Peak_Ground_Ratio_Upper+Peak_Ground_Ratio_Lower)/2;
				
		
				//Average value of the NPS Ground
				//NPS RADIAL PROFILE
				for (z=0; z<Ver_Profile.length; z++)
				{
					NPS_rad[z] = (rad_lower[z] + rad_upper[z])/2;
					print(NPS_rad[z]);
				}
				if(context == 1 && ImageType_array[k] == 1)//MAMMO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Radial_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(context == 2 && ImageType_array[k] == 2)//TOMO-PROJ
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Radial_TOMO_PROJ"+mm+"mm"+mAs+"mAs"+angle+"degrees");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3 && ImageType_array[k] == 1)//COMBO-MAMMO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Radial_COMBO_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}			
				if(context == 3 && ImageType_array[k] == 2 && dbg_val == 1)//COMBO-PROJ
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Radial_COMBO_PROJ"+mm+"mm"+mAs+"mAs"+angle+"degrees");
					selectWindow("Log");
					run("Close");
				}
				//NPS HORIZONTAL PROFILE
				for (z=0; z<Ver_Profile.length; z++)
				{
					NPS_hor[z] = (hor_lower[z] + hor_upper[z])/2;
					print(NPS_hor[z]);
				}
				if(context == 1 && ImageType_array[k] == 1)//MAMMO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Horizontal_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(context == 2 && ImageType_array[k] == 2 && dbg_val == 1)//TOMO-PROJ
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Horizontal_PROJECTION"+mm+"mm"+mAs+"mAs"+angle+"degrees");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3 && ImageType_array[k] == 1)//COMBO-MAMMO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Horizontal_COMBO_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3 && ImageType_array[k] == 2 && dbg_val == 1)//COMBO-PROJ
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Horizontal_COMBO_PROJ"+mm+"mm"+mAs+"mAs"+angle+"degrees");
					selectWindow("Log");
					run("Close");
				}
				//NPS VERTICAL PROFILE
				for (z=0; z<Ver_Profile.length; z++)
				{
					NPS_ver[z] = (ver_lower[z] + ver_upper[z])/2;
					print(NPS_ver[z]);
				}
				if(context == 1 && ImageType_array[k] == 1)//MAMMO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Vertical_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(context == 2 && ImageType_array[k] == 2 && dbg_val == 1)//TOMO-PROJ
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Vertical_TOMO_PROJ"+mm+"mm"+mAs+"mAs"+angle+"degrees");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3 && ImageType_array[k] == 1)//COMBO-MAMMO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Vertical_COMBO_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3 && ImageType_array[k] == 2 && dbg_val == 1)//COMBO-PROJ
				{
					selectWindow("Log");
					saveAs("Text", myDir+"NPS_Vertical_COMBO_PROJ"+mm+"mm"+mAs+"mAs"+angle+"degrees");
					selectWindow("Log");
					run("Close");
				}
			}
			//LSF MAMMO or MAMMO COMBO
			if(LSF_cal ==1 && ImageType_array[k] ==1) //Only for MAMMO
			{
				//Detection of the tungsten wire
				LSF_Offset_X = 35/pxdim; //Position of the tungsten wire with respect to the Al object in mm on X direction
				LSF_Offset_Y = 10/pxdim; //Position of the tungsten wire with respect to the Al object in mm on Y direction
				Wire_detection_ROI = 200; //ROI dimension for tungsten wire detection 
				
				selectWindow(filename);
				makeRectangle(Xobject+LSF_Offset_X-0.5*Wire_detection_ROI, Yobject-LSF_Offset_Y-0.5*Wire_detection_ROI, Wire_detection_ROI, Wire_detection_ROI);
				roiManager("Add");
				roiManager("Select", roiManager("count")-1);
				roiManager("Rename", "ROI_LSF");
				run("Copy");
				newImage("copy_image", "16-bit Black", Wire_detection_ROI, Wire_detection_ROI, 1);
				run("Paste");
				run("Invert");
				run("Enhance Contrast", "saturated=0.35");
				run("Select All");
				getRawStatistics(area, mean, min, max, std);
				Mean_Wire_detection_ROI = mean;
				//run("Subtract...", "value=Mean_Wire_detection_ROI");
				run("Fast Filters", "  filter=[background from median] x=10 y=0 preprocessing=smooth ");
				run("Enhance Contrast", "saturated=0.35");

				//Detection of the peak and coordinates calculation
				x = newArray(Wire_detection_ROI);
				Peak = newArray(Wire_detection_ROI);
				y = newArray(Wire_detection_ROI); 
				for(i = 0; i < Wire_detection_ROI; i++)
				{
					x[i] = Wire_detection_ROI-i;
					makeLine(i,0,i,Wire_detection_ROI);
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
				if (dbg_val ==1 )
				{
					selectWindow("Log");
					saveAs("Text", myDir+"Coordinates wire Left-Right direction"+mm+"mm"+mAs+"mAs");
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
				
				//Perpendicular Profiles of the wire 
				selectWindow(filename);
				run("Enhance Contrast", "saturated=0.35");
				makeRectangle(Xobject+LSF_Offset_X-0.5*Wire_detection_ROI, Yobject-LSF_Offset_Y-0.5*Wire_detection_ROI, Wire_detection_ROI, Wire_detection_ROI);
				run("Copy");
				newImage("copy_image", "16-bit Black", Wire_detection_ROI, Wire_detection_ROI, 1);
				run("Paste");
				run("Enhance Contrast", "saturated=0.35");
				selectWindow("copy_image");
				run("Invert");
				LSF_length = 128; //length of the wired profile -1 (real length 129);
				a_array = newArray(Wire_detection_ROI-20); //array containing the value of the gaussian fit parameters a
				b_array = newArray(Wire_detection_ROI-20); //array containing the value of the gaussian fit parameters b
				mu_array = newArray(Wire_detection_ROI-20); //array containing the value of the gaussian fit parameters mu
				sigma_array = newArray(Wire_detection_ROI-20); //array containing the value of the gaussian fit parameters sigma
				R_square_array = newArray(Wire_detection_ROI-20); //array containing the value of the gaussian fit parameters R_square
				lsf_left_right_sum = newArray(LSF_length);
				lsf_left_right_media = newArray(LSF_length);
				X_lsf = newArray(LSF_length);

				for(i = 10; i < Wire_detection_ROI-10; i++)
				{
					
					x[i] = i;
					XStart = x[i]+LSF_length/2*sin(angle_wire);
					YStart = y[i]-LSF_length/2*cos(angle_wire);
					XEnd = x[i]-LSF_length/2*sin(angle_wire);
					YEnd = y[i]+LSF_length/2*cos(angle_wire);
					run("Line Width...", "line=5");
					makeLine(XStart,YStart,XEnd,YEnd);
					roiManager("Add");
					roiManager("Select", roiManager("count")-1);
					roiManager("Rename", "ROI_LSF_profile"+i-9);
					lsf_left_right = getProfile();
					for (u=0; u < lsf_left_right.length; u++)
					{
						lsf_left_right_sum[u] += lsf_left_right[u];
					}
				
					for(j=0; j < lsf_left_right.length; j++)
					{
					
						X_lsf[j] = (j+1)*pxdim; //mm
					}
					
					Fit.doFit("gaussian", X_lsf,lsf_left_right);
					R_square = Fit.rSquared;
					a = Fit.p(0);
					b = Fit.p(1);
					mu = Fit.p(2);
					sigma = Fit.p(3);
					a_array[i-10] = a;
					b_array[i-10] = b;
					mu_array[i-10] = mu;
					sigma_array[i-10] = abs(sigma);
					R_square_array[i-10] = R_square;
					
					sigma_new = newArray(Wire_detection_ROI-20);
					a_new = newArray(Wire_detection_ROI-20);
					b_new = newArray(Wire_detection_ROI-20);
					mu_new = newArray(Wire_detection_ROI-20);
					new_count = 0;
					//Exclude the values for R^2 <0.3
					for(e=0; e < sigma_new.length; e++)
					{
						if(R_square_array[e] > 0.3)
						{
							a_new[new_count] = a_array[e];
							b_new[new_count] = b_array[e];
							mu_new[new_count] = mu_array[e];
							sigma_new[new_count] = sigma_array[e];
							new_count++;
						}
					}
				}
				//LSF average profile calculation and printing
				if(dbg_val ==1) print("mm\tPixel Value");
				for (u=0; u < lsf_left_right.length; u++)
				{
					lsf_left_right_media[u] = lsf_left_right_sum[u]/(Wire_detection_ROI-20);
					if(dbg_val ==1) print(d2s(X_lsf[u],4)+"\t"+d2s(lsf_left_right_media[u],2));
				}
				Fit.doFit("gaussian", X_lsf,lsf_left_right_media);
				R_square_gauss_mean = Fit.rSquared;
				a_gauss_mean = Fit.p(0);
				b_gauss_mean = Fit.p(1);
				mu_gauss_mean = Fit.p(2);
				sigma_gauss_mean = Fit.p(3);
				if(dbg_val ==1)
				{
					print("mm"+d2s(mm,2));
					print("mAs"+d2s(mAs,2));
					print("a"+d2s(a_gauss_mean,2));
					print("b"+d2s(b_gauss_mean,2));
					print("mu"+d2s(mu_gauss_mean,2));
					print("sigma"+d2s(sigma_gauss_mean,3));
					print("R^2"+d2s(R_square_gauss_mean,3));
				}
				if(context == 1 && dbg_val ==1)//MAMMO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"LSF_mean_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3 && dbg_val ==1)//COMBO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"LSF_mean_COMBO_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				a_stat = Array.trim(a_new, new_count);
				b_stat = Array.trim(b_new, new_count);
				mu_stat = Array.trim(mu_new, new_count);
				sigma_stat = Array.trim(sigma_new, new_count);
				
				Array.getStatistics(a_stat, min, max, mean, std);
				a_mean = mean;
				a_std = std;
				
				Array.getStatistics(b_stat, min, max, mean, std);
				b_mean = mean;
				b_std = std;
				
				Array.getStatistics(mu_stat, min, max, mean, std);
				mu_mean = mean;
				mu_std = std;

				Array.getStatistics(sigma_stat, min, max, mean, std);
				sigma_mean = mean;
				sigma_std = std;

				stat_data_number = sigma_stat.length;
				
				scale = 2*sqrt(2*log(2));
				lsf_array[k] = abs(sigma_gauss_mean*scale);
				if (dbg_val == 1)
				{
					print("Linear Fit parameters:\t");
					print("m:\t"+d2s(m,3));
					print("q:\t"+d2s(q,3));
					print("R2:\t"+d2s(R_square_lin,3));
					print("angle (degree):\t"+d2s(angle_wire_degree,3));
					print("Average value of gaussian fit parameters:\t");
					
					print("Used"+d2s(stat_data_number,0)+"profiles");
					
					print("a:\t"+d2s(a_mean,0)+"\t"+d2s(a_std,0));
					print("b:\t"+d2s(b_mean,-0)+"\t"+d2s(b_std,0));
					print("mu:\t"+d2s(mu_mean,3)+"\t"+d2s(mu_std,3));
					print("sigma mean:\t"+d2s(sigma_mean,3)+"\t"+d2s(sigma_std,3));
					
					print("Fit parameter for Average Gaussian curve:");
					print("a:\t"+d2s(a_gauss_mean,0));
					print("b:\t"+d2s(b_gauss_mean,-0));
					print("mu:\t"+d2s(mu_gauss_mean,3));
					print("sigma:\t"+d2s(sigma_gauss_mean,3));
					print("R2:\t"+d2s(R_square_gauss_mean,3));
				}
				if (context == 1 && dbg_val == 1)//MAMMO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"LSF_Result_Left-Right_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if (context == 3 && dbg_val == 1)//COMBO
				{
					selectWindow("Log");
					saveAs("Text", myDir+"LSF_Result_Left-Right_COMBO_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				
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
				makeRectangle(Xobject+LSF_Offset_X-0.5*Wire_detection_ROI, Yobject-LSF_Offset_Y-0.5*Wire_detection_ROI, Wire_detection_ROI, Wire_detection_ROI);
				run("Enhance Contrast", "saturated=0.35");
				roiManager("Add");
				roiManager("Select", roiManager("count")-1); 
				roiManager("Rename", "Pre_samp_LSF_ROI");
				run("Copy");
				newImage("Pre_sampling_LSF_ROI", "16-bit Black", Wire_detection_ROI, Wire_detection_ROI, 1);
				run("Paste");
				run("Invert");
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

					if(dbg_val == 1 && context == 1)
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
						saveAs("Text", myDir+"Pre_sampling_LSF_raw_MAMMO"+z+1+"lsf"+mm+"mm"+mAs+"mAs");
						selectWindow("Log");
						run("Close");
					}
					if(dbg_val == 1 && context == 3)
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
						saveAs("Text", myDir+"Pre_sampling_LSF_raw_COMBO_MAMMO"+z+1+"lsf"+mm+"mm"+mAs+"mAs");
						selectWindow("Log");
						run("Close");
					}					
					LSF_OC_stat_a[z] = a;
					LSF_OC_stat_b[z] = b;
					LSF_OC_stat_c[z] = c;
					LSF_OC_stat_d[z] = abs(d);
					LSF_OC_stat_FHWM[z] = abs(scale*d);					//vector containing the sigma value of the pre-sampling lsf
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

					//0 value for the pixel under 20 % of the peak and gaussian fit value for the tails
					for(j=0; j<lsfOC.length; j++)
					{
						if(lsfOC[j] <= 0.2)
						{
							lsfOC[j] = 0; //Put at 0 the ground value
						}
						if(lsfOC[j] > 0.2 && lsfOC[j] < 0.3)
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
					if(xOC_sort.length >= 2048)
					{
						slice_number = 1024;
					}
					if(xOC_sort.length >= 1024 && xOC_sort.length <= 2048)
					{
						slice_number = 512;
					}
					if(xOC_sort.length >= 512 && xOC_sort.length <= 1024)
					{
						slice_number = 256;
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
					if(context == 1 && dbg_val == 1)
					{
						print("a\t"+a);
						print("mu\t"+b);
						print("sigma\t"+c);
						selectWindow("Log");
						saveAs("Text", myDir+"Pre_sampling_LSF_normalized_buffer_MAMMO"+z+1+"lsf"+mm+"mm"+mAs+"mAs");
						selectWindow("Log");
						run("Close");
					}
					if(context == 3 && dbg_val == 1)
					{
						print("a\t"+a);
						print("mu\t"+b);
						print("sigma\t"+c);
						selectWindow("Log");
						saveAs("Text", myDir+"Pre_sampling_LSF_normalized_buffer_COMBO_MAMMO"+z+1+"lsf"+mm+"mm"+mAs+"mAs");
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
					if(dbg_val == 1 && context == 1)
					{						
						selectWindow("Log");
						saveAs("Text", myDir+"MTF_MAMMO"+(z+1)+"_"+mm+"mm"+d2s(mAs,0)+"mAs");
						selectWindow("Log");
						run("Close");
					}
					if(dbg_val == 1 && context == 3)
					{						
						selectWindow("Log");
						saveAs("Text", myDir+"MTF_COMBO_MAMMO"+(z+1)+"_"+mm+"mm"+d2s(mAs,0)+"mAs");
						selectWindow("Log");
						run("Close");
					}					
				}
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
				if(context == 1 && dbg_val == 1) //MAMMO
				{
					for(i=0; i<lsf_buffer.length; i++)
					{
						print(d2s(x_to_print[i],3)+"\t"+d2s(lsf_to_print[i],2)+"\t"+d2s(a_lsf_to_print*exp(-((x_to_print[i]-b_lsf_to_print)*(x_to_print[i]-b_lsf_to_print))/(2*c_lsf_to_print*c_lsf_to_print)),2));
					}
					print("a\t"+d2s(a_lsf_to_print,2));
					print("b\t"+d2s(b_lsf_to_print,2));
					print("c\t"+d2s(c_lsf_to_print,2));
					selectWindow("Log");
					saveAs("Text", myDir+"LSF_MAMMO"+mm+"mm"+d2s(mAs,0)+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3 && dbg_val == 1) //COMBO MAMMO
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
					saveAs("Text", myDir+"LSF_COMBO_MAMMO"+mm+"mm"+d2s(mAs,0)+"mAs");
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
				if(context == 1)
				{					
					selectWindow("Log");
					saveAs("Text", myDir+"MTF_MAMMO"+mm+"mm"+d2s(mAs,0)+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3)
				{					
					selectWindow("Log");
					saveAs("Text", myDir+"MTF_COMBO_MAMMO"+mm+"mm"+d2s(mAs,0)+"mAs");
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
				if(context == 1)
				{
					selectWindow("Log");
					saveAs("Text", myDir+"Pre_sampling_LSF_Fit_Results_MAMMO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
				if(context == 3)
				{
					selectWindow("Log");
					saveAs("Text", myDir+"Pre_sampling_LSF_Fit_Results_COMBO"+mm+"mm"+mAs+"mAs");
					selectWindow("Log");
					run("Close");
				}
			}
				
			//Test Results
			//Exposure Time Result
			if (ImageType_array[k] == 1 && exp_time_array[k] < 2000)//MAMMO
			{
				result_exp_time = "Pass";
			}
			if (ImageType_array[k] == 2 && exp_time_array[k] < 500)//PROJECTION
			{
				result_exp_time = "Pass";
			}
			//SNR Result
			if (ImageType_array[k] == 1)//MAMMO
			{
				if (Anode==" TUNGSTEN" && Filter==" SILVER" && SNR_BKG > 60)
				{
					result_SNR = "Pass";
				}
				if (Anode==" TUNGSTEN" && Filter==" RHODIUM " && SNR_BKG > 50)
				{
					result_SNR = "Pass";
				}
				if (Anode==" MOLYBDENUM" && Filter==" MOLYBDENUM" && SNR_BKG > 40)
				{
					result_SNR = "Pass";
				}
				if (Anode==" MOLYBDENUM" && Filter==" RHODIUM " && SNR_BKG > 40)
				{
					result_SNR = "Pass";
				}
				result_SNR_array[k] = result_SNR;
			}
			
			if (ImageType_array[k] == 2)//PROJECTION
			{
				if (Anode==" TUNGSTEN" && Filter==" SILVER" && SNR_BKG > 25)
				{
					result_SNR = "Pass";
				}
				if (Anode==" TUNGSTEN" && Filter==" RHODIUM " && SNR_BKG > 25)
				{
					result_SNR = "Pass";
				}
				result_SNR_array[k] = result_SNR;
			}

			//Homogeneity Result and ROISet saving
			if (ImageType_array[k] == 1)//MAMMO
			{
				if ((mm >35) && (mm <45))
				{
					if (Max_Mean_Dev[k] < 0.05)
					{
						result_Homogeneity_Mean = "Pass";
					}
					if (Max_Mean_Dev[k] < 0.15)
					{
						result_Homogeneity_SNR = "Pass";
					}
				}
			if(dbg_val ==1) roiManager("Save", myDir+"RoiSet_MAMMO"+mm+"mm"+mAs+"mAs"+".zip");
			roiManager("Delete");
			}

			if (ImageType_array[k] == 2)//PROJECTION
			{
				if ((mm >35) && (mm <45))
				{
					if (Max_Mean_Dev[k] < 0.1)
					{
						result_Homogeneity_Mean = "Pass";
					}
					if (Max_Mean_Dev[k] < 0.25)
					{
						result_Homogeneity_SNR = "Pass";
					}
				}
				if(context == 2)//TOMO
				{
				if(dbg_val ==1) roiManager("Save", myDir+"RoiSet_TOMO_PROJ"+mm+"mm"+mAs+"mAs"+angle+"degree"+".zip");
				//roiManager("Delete");
				}
				if(context == 3)//COMBO
				{
				if(dbg_val ==1) roiManager("Save", myDir+"RoiSet_COMBO_PROJ"+mm+"mm"+mAs+"mAs"+angle+"degree"+".zip");
				//roiManager("Delete");
				}
			}	

		}
	}


	run("Close All");
	roiManager("Reset");
	selectWindow("ROI Manager");
	run("Close");

}

//Bubble sort algorithm for the result printing (sort by thickness mm_array)
			scambio = 1;
			while(scambio)
			{
				scambio = 0;
				for(i = 0; i < list.length-1; i++)
				{
					if ((mm_array[i] > mm_array[i+1]) ||
						((mm_array[i] == mm_array[i+1]) && (angle_array[i] > angle_array[i+1])))
					{
						temp = mm_array[i];
						mm_array[i] = mm_array[i+1];
						mm_array[i+1] = temp;
						
						temp = angle_array[i];
						angle_array[i] = angle_array[i+1];
						angle_array[i+1] = temp;
						
						temp = Filter_array[i];
						Filter_array[i] = Filter_array[i+1];
						Filter_array[i+1] = temp;
						
						temp = kV_array[i];
						kV_array[i] = kV_array[i+1];
						kV_array[i+1] = temp;
						
						temp = mAs_array[i];
						mAs_array[i] = mAs_array[i+1];
						mAs_array[i+1] = temp;
						
						temp = uAs_array[i];
						uAs_array[i] = uAs_array[i+1];
						uAs_array[i+1] = temp;
						
						temp = mAs_tot_array[i];
						mAs_tot_array[i] = mAs_tot_array[i+1];
						mAs_tot_array[i+1] = temp;
						
						temp = exp_time_array[i];
						exp_time_array[i] = exp_time_array[i+1];
						exp_time_array[i+1] = temp;
						
						temp = mA_array[i];
						mA_array[i] = mA_array[i+1];
						mA_array[i+1] = temp;

						temp = ESAK_array[i];
						ESAK_array[i] = ESAK_array[i+1];
						ESAK_array[i+1] = temp;
	
						temp = HVL_array[i];
						HVL_array[i] = HVL_array[i+1];
						HVL_array[i+1] = temp;

						temp = MGD_array[i];
						MGD_array[i] = MGD_array[i+1];
						MGD_array[i+1] = temp;
						
						temp = ImageType_array[i];
						ImageType_array[i] = ImageType_array[i+1];
						ImageType_array[i+1] = temp;

						temp = MeanBKG_array[i];
						MeanBKG_array[i] = MeanBKG_array[i+1];
						MeanBKG_array[i+1] = temp;

						temp = StdBKG_array[i];
						StdBKG_array[i] = StdBKG_array[i+1];
						StdBKG_array[i+1] = temp;

						temp = SNR_BKG_array[i];
						SNR_BKG_array[i] = SNR_BKG_array[i+1];
						SNR_BKG_array[i+1] = temp;
						
						temp = MeanDetail_array[i];
						MeanDetail_array[i] = MeanDetail_array[i+1];
						MeanDetail_array[i+1] = temp;
						
						temp = StdDetail_array[i];
						StdDetail_array[i] = StdDetail_array[i+1];
						StdDetail_array[i+1] = temp;
						
						temp = SNR_Object_array[i];
						SNR_Object_array[i] = SNR_Object_array[i+1];
						SNR_Object_array[i+1] = temp;
						
						temp = result_SNR_array[i];
						result_SNR_array[i] = result_SNR_array[i+1];
						result_SNR_array[i+1] = temp;
						
						temp = Contrast_array[i];
						Contrast_array[i] = Contrast_array[i+1];
						Contrast_array[i+1] = temp;

						temp = CNR_array[i];
						CNR_array[i] = CNR_array[i+1];
						CNR_array[i+1] = temp;
						
						temp = MeanROI1[i];
						MeanROI1[i] = MeanROI1[i+1];
						MeanROI1[i+1] = temp;
						
						temp = MeanROI2[i];
						MeanROI2[i] = MeanROI2[i+1];
						MeanROI2[i+1] = temp;

						temp = MeanROI3[i];
						MeanROI3[i] = MeanROI3[i+1];
						MeanROI3[i+1] = temp;

						temp = MeanROI4[i];
						MeanROI4[i] = MeanROI4[i+1];
						MeanROI4[i+1] = temp;
						
						temp = MeanROI5[i];
						MeanROI5[i] = MeanROI5[i+1];
						MeanROI5[i+1] = temp;
						
						temp = SNR_ROI1[i];
						SNR_ROI1[i] = SNR_ROI1[i+1];
						SNR_ROI1[i+1] = temp;

						temp = SNR_ROI2[i];
						SNR_ROI2[i] = SNR_ROI2[i+1];
						SNR_ROI2[i+1] = temp;

						temp = SNR_ROI3[i];
						SNR_ROI3[i] = SNR_ROI3[i+1];
						SNR_ROI3[i+1] = temp;

						temp = SNR_ROI4[i];
						SNR_ROI4[i] = SNR_ROI4[i+1];
						SNR_ROI4[i+1] = temp;

						temp = SNR_ROI5[i];
						SNR_ROI5[i] = SNR_ROI5[i+1];
						SNR_ROI5[i+1] = temp;

						temp = Mean_Dev1[i];
						Mean_Dev1[i] = Mean_Dev1[i+1];
						Mean_Dev1[i+1] = temp;

						temp = Mean_Dev2[i];
						Mean_Dev2[i] = Mean_Dev2[i+1];
						Mean_Dev2[i+1] = temp;

						temp = Mean_Dev3[i];
						Mean_Dev3[i] = Mean_Dev3[i+1];
						Mean_Dev3[i+1] = temp;

						temp = Mean_Dev4[i];
						Mean_Dev4[i] = Mean_Dev4[i+1];
						Mean_Dev4[i+1] = temp;

						temp = Mean_Dev5[i];
						Mean_Dev5[i] = Mean_Dev5[i+1];
						Mean_Dev5[i+1] = temp;

						temp = SNR_Dev1[i];
						SNR_Dev1[i] = SNR_Dev1[i+1];
						SNR_Dev1[i+1] = temp;

						temp = SNR_Dev2[i];
						SNR_Dev2[i] = SNR_Dev2[i+1];
						SNR_Dev2[i+1] = temp;

						temp = SNR_Dev3[i];
						SNR_Dev3[i] = SNR_Dev3[i+1];
						SNR_Dev3[i+1] = temp;

						temp = SNR_Dev4[i];
						SNR_Dev4[i] = SNR_Dev4[i+1];
						SNR_Dev4[i+1] = temp;

						temp = SNR_Dev5[i];
						SNR_Dev5[i] = SNR_Dev5[i+1];
						SNR_Dev5[i+1] = temp;

						temp = MeanTot[i];
						MeanTot[i] = MeanTot[i+1];
						MeanTot[i+1] = temp;

						temp = SNRTot[i];
						SNRTot[i] = SNRTot[i+1];
						SNRTot[i+1] = temp;

						temp = Max_Mean_Dev[i];
						Max_Mean_Dev[i] = Max_Mean_Dev[i+1];
						Max_Mean_Dev[i+1] = temp;

						temp = Max_SNR_Dev[i];
						Max_SNR_Dev[i] = Max_SNR_Dev[i+1];
						Max_SNR_Dev[i+1] = temp;

						temp = grid_visibility_array[i];
						grid_visibility_array[i] = grid_visibility_array[i+1];
						grid_visibility_array[i+1] = temp;

						temp = Grid_Peak_array[i];
						Grid_Peak_array[i] = Grid_Peak_array[i+1];
						Grid_Peak_array[i+1] = temp;

						temp = NPS_Ground_array[i];
						NPS_Ground_array[i] = NPS_Ground_array[i+1];
						NPS_Ground_array[i+1] = temp;

						temp = BackGround_mean_array[i];
						BackGround_mean_array[i] = BackGround_mean_array[i+1];
						BackGround_mean_array[i+1] = temp;
						scambio = 1;
						
						temp = BackGround_std_array[i];
						BackGround_std_array[i] = BackGround_std_array[i+1];
						BackGround_std_array[i+1] = temp;

						temp = BackGround_SNR_array[i];
						BackGround_SNR_array[i] = BackGround_SNR_array[i+1];
						BackGround_SNR_array[i+1] = temp;

						temp = Al_mean_array[i];
						Al_mean_array[i] = Al_mean_array[i+1];
						Al_mean_array[i+1] = temp;
						
						temp = Al_std_array[i];
						Al_std_array[i] = Al_std_array[i+1];
						Al_std_array[i+1] = temp;

						temp = Al_SNR_array[i];
						Al_SNR_array[i] = Al_SNR_array[i+1];
						Al_SNR_array[i+1] = temp;

						temp = Contrast_rec_array[i];
						Contrast_rec_array[i] = Contrast_rec_array[i+1];
						Contrast_rec_array[i+1] = temp;

						temp = CNR_rec_array[i];
						CNR_rec_array[i] = CNR_rec_array[i+1];
						CNR_rec_array[i+1] = temp;

						temp = sigma_asf[i];
						sigma_asf[i] = sigma_asf[i+1];
						sigma_asf[i+1] = temp;

						temp = CNR_asf_focused_slice[i];
						CNR_asf_focused_slice[i] = CNR_asf_focused_slice[i+1];
						CNR_asf_focused_slice[i+1] = temp;
						
						temp = lsf_array[i];
						lsf_array[i] = lsf_array[i+1];
						lsf_array[i+1] = temp;
						
						temp = wire_angle_array[i];
						wire_angle_array[i] = wire_angle_array[i+1];
						wire_angle_array[i+1] = temp;

						temp = chest_nipple_alingment_results[i];
						chest_nipple_alingment_results[i] = chest_nipple_alingment_results[i+1];
						chest_nipple_alingment_results[i+1] = temp;

						temp = left_right_alingment_results[i];
						left_right_alingment_results[i] = left_right_alingment_results[i+1];
						left_right_alingment_results[i+1] = temp;

						temp = focused_slice_chest_array[i];
						focused_slice_chest_array[i] = focused_slice_chest_array[i+1];
						focused_slice_chest_array[i+1] = temp;

						temp = focused_slice_nipple_array[i];
						focused_slice_nipple_array[i] = focused_slice_nipple_array[i+1];
						focused_slice_nipple_array[i+1] = temp;

						temp = focused_slice_left_array[i];
						focused_slice_left_array[i] = focused_slice_left_array[i+1];
						focused_slice_left_array[i+1] = temp;

						temp = focused_slice_right_array[i];
						focused_slice_right_array[i] = focused_slice_right_array[i+1];
						focused_slice_right_array[i+1] = temp;

						temp = SliceNumber_array[i];
						SliceNumber_array[i] = SliceNumber_array[i+1];
						SliceNumber_array[i+1] = temp;
						scambio = 1;
					}
				}
			}
			
//Print of the output file for the AEC Test
if (context == 1)//MAMMO
{
	print("MAMMO CONTEXT: Automatic Exposure Control (AEC) system Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("\nmm\tMODE\tANODE\tFILTER\tkV\tmAs\tExp. Time (ms)\tAnodic Current (mA)\tESAK (mGy)\tHVL (mmAl)\tMGD (mGy)\tDose Rate (mGy/s)");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1)//MAMMO
		{
			print(d2s(mm_array[j],0),"\t"+Exp_mode+"\t"+Anode+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(exp_time_array[j],0)+"\t"+d2s(mA_array[j],0)+"\t"+d2s(ESAK_array[j],2)+"\t"+d2s(HVL_array[j],3)+"\t"+d2s(MGD_array[j]*100,2)+"\t"+d2s((ESAK_array[j])/((exp_time_array[j])/1000),2));
		}
	}
	if (result_exp_time == "Pass")
	{
		print("Tube Efficiency OK. Exposure time for 40 mm is under the limit of 1s");
	}
	if (result_exp_time == "Fail")
	{
		print("Tube Efficiency Fail. Exposure time for 40 mm is over the limit of 1s");
	}
	selectWindow("Log");
	saveAs("Text", myDir+"AEC_MAMMO.txt");
	selectWindow("Log");
	run("Close");
}
if (context == 2)//TOMO
{
	print("TOMO CONTEXT PROJECTION images: Automatic Exposure Control (AEC) system Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("\nmm\tAngle (°)\tMODE\tANODE\tFILTER\tkV\tmAs\tExp. Time (ms)\tAnodic Current (mA)\tESAK (mGy)\tHVL (mmAl)\tMGD (mGy)\tDose Rate (mGy/s)");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 2)//PROJECTION
		{
			print(d2s(mm_array[j],0)+"\t"+d2s(angle_array[j],2)+"\t"+Exp_mode+"\t"+Anode+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(uAs_array[j]/1000,1)+"\t"+d2s(exp_time_array[j],0)+"\t"+d2s(mA_array[j],0)+"\t"+d2s(ESAK_array[j],2)+"\t"+d2s(HVL_array[j],3)+"\t"+d2s(MGD_array[j]*100,2)+"\t"+d2s((ESAK_array[j])/((exp_time_array[j])/1000),2));
		}		
	}
	if (result_exp_time == "Pass")
	{
		print("Tube Efficiency OK. Exposure time for 40 mm is under the limit of 1s");
	}
	if (result_exp_time == "Fail")
	{
		print("Tube Efficiency Fail. Exposure time for 40 mm is over the limit of 1s");
	}
	selectWindow("Log");
	saveAs("Text", myDir+"AEC_TOMO_Projection.txt");
	selectWindow("Log");
	run("Close");
}
if (context == 3)//COMBO
{
	print("COMBO CONTEXT MAMMO images: Automatic Exposure Control (AEC) system Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("\nmm\tMODE\tANODE\tFILTER\tkV\tmAs\tExp. Time (ms)\tAnodic Current (mA)\tESAK (mGy)\tHVL (mmAl)\tMGD (mGy)\tDose Rate (mGy/s)");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1)//MAMMO
		{
			print(d2s(mm_array[j],0),"\t"+Exp_mode+"\t"+Anode+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(exp_time_array[j],0)+"\t"+d2s(mA_array[j],0)+"\t"+d2s(ESAK_array[j],2)+"\t"+d2s(HVL_array[j],3)+"\t"+d2s(MGD_array[j]*100,2)+"\t"+d2s((ESAK_array[j])/((exp_time_array[j])/1000),2));
		}
	}
	if (result_exp_time == "Pass")
	{
		print("Tube Efficiency OK. Exposure time for 40 mm is under the limit of 1s");
	}
	if (result_exp_time == "Fail")
	{
		print("Tube Efficiency Fail. Exposure time for 40 mm is over the limit of 1s");
	}
	selectWindow("Log");
	saveAs("Text", myDir+"AEC_MAMMO_COMBO.txt");
	selectWindow("Log");
	run("Close");
	
	print("COMBO CONTEXT PROJECTION images: Automatic Exposure Control (AEC) system Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("\nmm\tAngle (°)\tMODE\tANODE\tFILTER\tkV\tmAs\tExp. Time (ms)\tAnodic Current (mA)\tESAK (mGy)\tHVL (mmAl)\tMGD (mGy)\tDose Rate (mGy/s)");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 2)//PROJECTION
		{
			print(d2s(mm_array[j],0)+"\t"+d2s(angle_array[j],2)+"\t"+Exp_mode+"\t"+Anode+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(uAs_array[j]/1000,2)+"\t"+d2s(exp_time_array[j],0)+"\t"+d2s(mA_array[j],0)+"\t"+d2s(ESAK_array[j],2)+"\t"+d2s(HVL_array[j],3)+"\t"+d2s(MGD_array[j]*100,2)+"\t"+d2s((ESAK_array[j])/((exp_time_array[j])/1000),2));
		}		
	}
	if (result_exp_time == "Pass")
	{
		print("Tube Efficiency OK. Exposure time for 40 mm is under the limit of 1s");
	}
	if (result_exp_time == "Fail")
	{
		print("Tube Efficiency Fail. Exposure time for 40 mm is over the limit of 1s");
	}
	selectWindow("Log");
	saveAs("Text", myDir+"AEC_COMBO_Projection.txt");
	selectWindow("Log");
	run("Close");
}
//Print of the output file for the SNR-CNR Test
if (context == 1)//MAMMO
{
	print("MAMMO CONTEXT: Signal-to-Noise Ration (SNR) and Contrast-to-Noise Ratio (CNR) Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("\nmm\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tLSF\tSNR Result");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1)//MAMMO
		{
		   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(MeanBKG_array[j],1)+"\t"+d2s(StdBKG_array[j],2)+"\t"+d2s(SNR_BKG_array[j],1)+"\t"+d2s(MeanDetail_array[j],1)+"\t"+d2s(StdDetail_array[j],2)+"\t"+d2s(SNR_Object_array[j],1)+"\t"+d2s(Contrast_array[j],3)+"\t"+d2s(CNR_array[j],1)+"\t"+d2s(lsf_array[j],3)+"\t"+result_SNR_array[j]);
		}
	}
	print("\n\n\n\n\n\n\nTest Limit SNR > 60");
	selectWindow("Log");
	saveAs("Text", myDir+"CNR_MAMMO.txt");
	selectWindow("Log");
	run("Close");
}
if (context == 2)//TOMO
{
	print("TOMO CONTEXT PROJECTION images: Signal-to-Noise Ration (SNR) and Contrast-to-Noise Ratio (CNR) Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("mm\tAngle (°)\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tSNR Result\tCNR ASF");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 2)//PROJECTION
		{
			print(d2s(mm_array[j],0)+"\t"+d2s(angle_array[j],2)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(uAs_array[j]/1000,2)+"\t"+d2s(MeanBKG_array[j],1)+"\t"+d2s(StdBKG_array[j],2)+"\t"+d2s(SNR_BKG_array[j],1)+"\t"+d2s(MeanDetail_array[j],1)+"\t"+d2s(StdDetail_array[j],2)+"\t"+d2s(SNR_Object_array[j],1)+"\t"+d2s(Contrast_array[j],3)+"\t"+d2s(CNR_array[j],2)+"\t"+result_SNR_array[j]+"\t"+d2s(CNR_asf_focused_slice[j],2));
		}
	}
	print("\nTest Limit SNR > 25");
	selectWindow("Log");
	saveAs("Text", myDir+"CNR_TOMO_Projection.txt");
	selectWindow("Log");
	run("Close");
	
	print("TOMO CONTEXT RECONSTRUCTED images: Signal-to-Noise Ration (SNR) and Contrast-to-Noise Ratio (CNR) Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("mm\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tsigma LSF\tSNR Result\tCNR ASF");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3)//RECO
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(BackGround_mean_array[j],1)+"\t"+d2s(BackGround_std_array[j],2)+"\t"+d2s(BackGround_SNR_array[j],1)+"\t"+d2s(Al_mean_array[j],1)+"\t"+d2s(Al_std_array[j],2)+"\t"+d2s(Al_SNR_array[j],1)+"\t"+d2s(Contrast_rec_array[j],3)+"\t"+d2s(CNR_asf_focused_slice[j],2)+"\t"+d2s(lsf_array[j],3)+"\t"+result_SNR_array[j]+"\t"+d2s(CNR_rec_array[j],2));
		}
	}
	print("\nTest Limit SNR > 50");
	selectWindow("Log");
	saveAs("Text", myDir+"CNR_TOMO_RECONSTRUCTION.txt");
	selectWindow("Log");
	run("Close");
}
if (context == 3)//COMBO
{
	print("COMBO CONTEXT MAMMO images: Signal-to-Noise Ration (SNR) and Contrast-to-Noise Ratio (CNR) Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("mm\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tLSF\tSNR Result");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1)//MAMMO
		{
		   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(MeanBKG_array[j],1)+"\t"+d2s(StdBKG_array[j],2)+"\t"+d2s(SNR_BKG_array[j],1)+"\t"+d2s(MeanDetail_array[j],1)+"\t"+d2s(StdDetail_array[j],2)+"\t"+d2s(SNR_Object_array[j],1)+"\t"+d2s(Contrast_array[j],3)+"\t"+d2s(CNR_array[j],1)+"\t"+d2s(lsf_array[j],3)+"\t"+result_SNR_array[j]);
		}
	}
	print("\nTest Limit SNR > 60");
	selectWindow("Log");
	saveAs("Text", myDir+"CNR_MAMMO_COMBO.txt");
	selectWindow("Log");
	run("Close");
	
	print("COMBO CONTEXT PROJECTION images:Signal-to-Noise Ration (SNR) and Contrast-to-Noise Ratio (CNR) Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(Distance_Source_to_patient,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("mm\tAngle (°)\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tSNR Result");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 2)//PROJECTION
		{
			print(d2s(mm_array[j],0)+"\t"+d2s(angle_array[j],2)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(uAs_array[j]/1000,1)+"\t"+d2s(MeanBKG_array[j],1)+"\t"+d2s(StdBKG_array[j],2)+"\t"+d2s(SNR_BKG_array[j],1)+"\t"+d2s(MeanDetail_array[j],1)+"\t"+d2s(StdDetail_array[j],2)+"\t"+d2s(SNR_Object_array[j],1)+"\t"+d2s(Contrast_array[j],3)+"\t"+d2s(CNR_array[j],1)+"\t"+result_SNR_array[j]);
		}
	}
	print("\nTest Limit SNR > 25");
	selectWindow("Log");
	saveAs("Text", myDir+"CNR_COMBO_Projection.txt");
	selectWindow("Log");
	run("Close");
	
	print("COMBO CONTEXT RECONSTRUCTED images: Signal-to-Noise Ration (SNR) and Contrast-to-Noise Ratio (CNR) Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("\nmmA\tFILTER\tkV\tmAs tot.\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tLSF\tSNR Result\tCNR ASF");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3)//RECO
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(BackGround_mean_array[j],1)+"\t"+d2s(BackGround_std_array[j],2)+"\t"+d2s(BackGround_SNR_array[j],1)+"\t"+d2s(Al_mean_array[j],1)+"\t"+d2s(Al_std_array[j],2)+"\t"+d2s(Al_SNR_array[j],1)+"\t"+d2s(Contrast_rec_array[j],3)+"\t"+d2s(CNR_asf_focused_slice[j],2)+"\t"+d2s(lsf_array[j],3)+"\t"+result_SNR_array[j]+"\t"+d2s(CNR_rec_array[j],2));
		}
	}
	print("\nTest Limit SNR > 50");
	selectWindow("Log");
	saveAs("Text", myDir+"CNR_COMBO_RECONSTRUCTION.txt");
	selectWindow("Log");
	run("Close");
}
//Print of the HOMOGENEITY Test result
if (context == 1)//MAMMO
{
	print("MAMMO CONTEXT: Homogeneity Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("\nmm\tFILTER\tkV\tmAs\tROI1\tROI2\tROI3\tROI4\tROI5\tAverage\tDev max % Mean");
	
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1) //MAMMO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(MeanROI1[j],1)+"\t"+d2s(MeanROI2[j],2)+"\t"+d2s(MeanROI3[j],1)+"\t"+d2s(MeanROI4[j],1)+"\t"+d2s(MeanROI5[j],2)+"\t"+d2s(MeanTot[j],1)+"\t"+d2s(Max_Mean_Dev[j],3));
		}
	}
	print("\nmm\tFILTER\tkV\tmAs\tSNR1\tSNR2\tSNR3\tSNR4\tSNR5\tAverage SNR\tDev % SNR");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1) //MAMMO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(SNR_ROI1[j],1)+"\t"+d2s(SNR_ROI2[j],1)+"\t"+d2s(SNR_ROI3[j],1)+"\t"+d2s(SNR_ROI4[j],1)+"\t"+d2s(SNR_ROI5[j],1)+"\t"+d2s(SNRTot[j],1)+"\t"+d2s(Max_SNR_Dev[j],3));
		}
	}
		print("\nHomogeneity Test Result for Mean Grey Level\t"+result_Homogeneity_Mean);
		print("Limit Dev % max for average grey level < 5%");
		print("Homogeneity Test Result for SNR\t"+result_Homogeneity_SNR);
		print("Limit Dev % max for SNR < 15%");
		selectWindow("Log");
		saveAs("Text", myDir+"Homogeneity_MAMMO.txt");
		selectWindow("Log");
		run("Close");
}
if (context == 2)//TOMO
{
	print("TOMO CONTEXT RECONSTRUCTED images: Homogeneity Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("\nmm\tFILTER\tkV\tmAs tot.\tROI1\tROI2\tROI3\tROI4\tROI5\tAverage\tDev max % Mean");
	
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3) //RECO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(MeanROI1[j],1)+"\t"+d2s(MeanROI2[j],2)+"\t"+d2s(MeanROI3[j],1)+"\t"+d2s(MeanROI4[j],1)+"\t"+d2s(MeanROI5[j],2)+"\t"+d2s(MeanTot[j],1)+"\t"+d2s(Max_Mean_Dev[j],3));
		}
	}
	print("\nmm\tFILTER\tkV\tmAs\tSNR1\tSNR2\tSNR3\tSNR4\tSNR5\tAverage SNR\tDev % SNR");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3) //RECO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(SNR_ROI1[j],1)+"\t"+d2s(SNR_ROI2[j],1)+"\t"+d2s(SNR_ROI3[j],1)+"\t"+d2s(SNR_ROI4[j],1)+"\t"+d2s(SNR_ROI5[j],1)+"\t"+d2s(SNRTot[j],1)+"\t"+d2s(Max_SNR_Dev[j],3));
		}
	}
		print("\nHomogeneity Test Result for Mean Grey Level\t"+result_Homogeneity_Mean);
		print("Limit Dev % max for average grey level < 5%");
		print("Homogeneity Test Result for SNR\t"+result_Homogeneity_SNR);
		print("Limit Dev % max for SNR < 15%");
		selectWindow("Log");
		saveAs("Text", myDir+"Homogeneity_TOMO_RECO.txt");
		selectWindow("Log");
		run("Close");
}
if (context == 3)//COMBO
{
	print("COMBO_MAMMO CONTEXT: Homogeneity Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("\nmm\tFILTER\tkV\tmAs\tROI1\tROI2\tROI3\tROI4\tROI5\tAverage\tDev max % Mean");
	
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1) //MAMMO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(MeanROI1[j],1)+"\t"+d2s(MeanROI2[j],2)+"\t"+d2s(MeanROI3[j],1)+"\t"+d2s(MeanROI4[j],1)+"\t"+d2s(MeanROI5[j],2)+"\t"+d2s(MeanTot[j],1)+"\t"+d2s(Max_Mean_Dev[j],3));
		}
	}
	print("\nmm\tFILTER\tkV\tmAs\tSNR1\tSNR2\tSNR3\tSNR4\tSNR5\tAverage SNR\tDev % SNR");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1) //MAMMO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(SNR_ROI1[j],1)+"\t"+d2s(SNR_ROI2[j],1)+"\t"+d2s(SNR_ROI3[j],1)+"\t"+d2s(SNR_ROI4[j],1)+"\t"+d2s(SNR_ROI5[j],1)+"\t"+d2s(SNRTot[j],1)+"\t"+d2s(Max_SNR_Dev[j],3));
		}
	}
	print("\nHomogeneity Test Result for Mean Grey Level\t"+result_Homogeneity_Mean);
	print("Limit Dev % max for average grey level < 5%");
	print("Homogeneity Test Result for SNR\t"+result_Homogeneity_SNR);
	print("Limit Dev % max for SNR < 15%");
	selectWindow("Log");
	saveAs("Text", myDir+"Homogeneity_COMBO_MAMMO.txt");
	selectWindow("Log");
	run("Close");
		
	print("COMBO CONTEXT RECONSTRUCTED images: Homogeneity Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("\nmm\tFILTER\tkV\tmAs tot.\tROI1\tROI2\tROI3\tROI4\tROI5\tAverage\tDev max % Mean");
	
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3) //RECO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(MeanROI1[j],1)+"\t"+d2s(MeanROI2[j],2)+"\t"+d2s(MeanROI3[j],1)+"\t"+d2s(MeanROI4[j],1)+"\t"+d2s(MeanROI5[j],2)+"\t"+d2s(MeanTot[j],1)+"\t"+d2s(Max_Mean_Dev[j],3));
		}
	}
	print("\n mm PMMA\t FILTER\t kV\t mAs\t SNR1\t SNR2\t SNR3\t SNR4\t SNR5\t Average SNR\t Dev % SNR");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3) //RECO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(SNR_ROI1[j],1)+"\t"+d2s(SNR_ROI2[j],1)+"\t"+d2s(SNR_ROI3[j],1)+"\t"+d2s(SNR_ROI4[j],1)+"\t"+d2s(SNR_ROI5[j],1)+"\t"+d2s(SNRTot[j],1)+"\t"+d2s(Max_SNR_Dev[j],3));
		}
	}
		print("\nHomogeneity Test Result for Mean Grey Level\t"+result_Homogeneity_Mean);
		print("Limit Dev % max for average grey level < 5%");
		print("Homogeneity Test Result for SNR\t"+result_Homogeneity_SNR);
		print("Limit Dev % max for SNR < 15%");
		selectWindow("Log");
		saveAs("Text", myDir+"Homogeneity_COMBO_RECO.txt");
		selectWindow("Log");
		run("Close");
}
//Print of the first line of output file with final results for GRID VISIBILITY-NOISE Test
if (context == 1)//MAMMO
{
	print("MAMMO CONTEXT: Grid Visibility-Noise Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("\nmm\tFILTER\tkV\tmAs\tExp. Time (ms)\tGrid Peak\tNPS Ground\tPeak-Ground Ratio");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1)//MAMMO
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(exp_time_array[j],0)+"\t"+d2s(Grid_Peak_array[j],-1)+"\t"+d2s(NPS_Ground_array[j],-1)+"\t"+d2s(grid_visibility_array[j],1));
		}
	}
	print("\nLimit Peak-Ground Ratio < 20");
	selectWindow("Log");
	saveAs("Text", myDir+"Grid_Visibility.txt");
	selectWindow("Log");
	run("Close");
}
if (context == 2)//TOMO
{
	print("TOMO CONTEXT: Projection Noise Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("\nmm\tangle\tFILTER\tkV\tmAs\tExp. Time (ms)\tNPS Ground");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 2)//PROJECTION
		{
		   print(d2s(mm_array[j],0)+"\t"+d2s(angle_array[j],2)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s((uAs_array[j])/1000,1)+"\t"+d2s(exp_time_array[j],0)+"\t"+d2s(NPS_Ground_array[j],-1));
		}
	}
	selectWindow("Log");
	saveAs("Text", myDir+"Noise_TOMO_Projection.txt");
	selectWindow("Log");
	run("Close");
	
	print("TOMO RECO CONTEXT: Noise Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("\nmm\tFILTER\tkV\tmAs tot.\tExp. Time (ms)\tNPS Ground");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3)//RECO
		{
		   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(NPS_Ground_array[j],-1));
		}
	}
	selectWindow("Log");
	saveAs("Text", myDir+"Noise_TOMO_RECO.txt");
	selectWindow("Log");
	run("Close");
}
if (context == 3)//COMBO
{
	print("COMBO CONTEXT PROJECTION images: Noise Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("\nmm\tFILTER\tkV\tmAs\tExp. Time (ms)\tNPS Ground");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 2)//PROJECTION
		{
		   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s((uAs_array[j])/1000,1)+"\t"+d2s(exp_time_array[j],0)+"\t"+d2s(NPS_Ground_array[j],-1));
		}
	}
	selectWindow("Log");
	saveAs("Text", myDir+"Noise_COMBO_Projection.txt");
	selectWindow("Log");
	run("Close");
	
	print("COMBO CONTEXT MAMMO images: Noise Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("\nmm\tFILTER\tkV\tmAs\tNPS Ground");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1)//MAMMO
		{
		   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(NPS_Ground_array[j],-1));
		   
		}
	}
	selectWindow("Log");
	saveAs("Text", myDir+"Noise_COMBO_MAMMO.txt");
	selectWindow("Log");
	run("Close");
	
	print("COMBO RECO CONTEXT: Noise Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	print("\nmm\tFILTER\tkV\tmAs tot.\tNPS Ground");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3)//RECO
		{
		   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(NPS_Ground_array[j],-1));
		}
	}
	selectWindow("Log");
	saveAs("Text", myDir+"Noise_COMBO_RECO.txt");
	selectWindow("Log");
	run("Close");
}
//Printing of the Alignment Test Results
if(context == 2 && LSF_cal == 1) //TOMO
{
	print("TOMO RECO CONTEXT: Alignment Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(Distance_Source_to_patient,1));
	print("\nmm\tFILTER\tkV\tmAs tot.\tFocused Slice at Chest side\tFocused Slice at Nipple side\tAlignment chest-nipple result\tFocused Slice at Right side\t Focused Slice at Left side\t Alignment left-right result\t");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3)//RECO
		{
		   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(focused_slice_chest_array[j],0)+"\t"+d2s(focused_slice_nipple_array[j],0)+"\t"+chest_nipple_alingment_results[j]+"\t"+d2s(focused_slice_right_array[j],0)+"\t"+d2s(focused_slice_left_array[j],0)+"\t"+d2s(SliceNumber_array[j],0)+"\t"+left_right_alingment_results[j]);
		}
	}
	selectWindow("Log");
	saveAs("Text", myDir+"Alignment_Test_TOMO.txt");
	selectWindow("Log");
	run("Close");
}
if(context == 3 && LSF_cal == 1) //COMBO
{
	print("TOMO COMBO CONTEXT: Alignment Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(Distance_Source_to_patient,1));
	print("\nmm\tFILTER\tkV\tmAs tot. \tFocused Slice at Chest side\tFocused Slice at Nipple side\tAlignment chest-nipple result\tFocused Slice at Right side\t Focused Slice at Left side\t Alignment left-right result\t");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3)//RECO
		{
		   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(focused_slice_chest_array[j],0)+"\t"+d2s(focused_slice_nipple_array[j],0)+"\t"+chest_nipple_alingment_results[j]+"\t"+d2s(focused_slice_right_array[j],0)+"\t"+d2s(focused_slice_left_array[j],0)+"\t"+d2s(SliceNumber_array[j],0)+"\t"+left_right_alingment_results[j]);
		}
	}
	selectWindow("Log");
	saveAs("Text", myDir+"Alignment_Test_COMBO.txt");
	selectWindow("Log");
	run("Close");
}
//Printing of the LSF-MTF Test Results
if(context == 1 && LSF_cal == 1) //MAMMO
{
	print("MAMMO CONTEXT: LSF Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("\nmm\tFILTER\tkV\tmAs\tangle (degree)\tFHWM (mm)");
	for(j = 0; j < list.length; j++)
	{
	   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(wire_angle_array[j],3)+"\t"+d2s(lsf_array[j],3));
	}
	selectWindow("Log");
	saveAs("Text", myDir+"LSF_Result_MAMMO.txt");
	selectWindow("Log");
	run("Close");
}
if(context == 2 && LSF_cal == 1) //TOMO
{
	print("TOMO CONTEXT: LSF Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("\nmmPMMA\tFILTER\tkV\tmAstot.\tangle (degree)\tFHWM (mm)");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3)//RECO
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(wire_angle_array[j],3)+"\t"+d2s(lsf_array[j],3));
		}
	}
	selectWindow("Log");
	saveAs("Text", myDir+"LSF_Result_TOMO.txt");
	selectWindow("Log");
	run("Close");
}
if(context == 3 && LSF_cal == 1) //COMBO MAMMO
{
	print("COMBO MAMMO CONTEXT: LSF Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("\nmmPMMA\tFILTER\tkV\tmAs\tangle (degree)\tFHWM (mm)");
	
	for(j = 0; j < list.length; j++)
	{
	   if (ImageType_array[j] == 1)//MAMMO
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(wire_angle_array[j],3)+"\t"+d2s(lsf_array[j],3));
		}
	}
	
	selectWindow("Log");
	saveAs("Text", myDir+"LSF_Result_COMBO_MAMMO.txt");
	selectWindow("Log");
	run("Close");
}
if(context == 3 && LSF_cal == 1 ) //COMBO RECO
{
	print("COMBO RECO CONTEXT: LSF Test");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("\nmmPMMA\tFILTER\tkV\tmAstot.\tangle (degree)\tFHWM (mm)");
	for(j = 0; j < list.length; j++)
	{
	   if (ImageType_array[j] == 3)//RECO
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(wire_angle_array[j],3)+"\t"+d2s(lsf_array[j],3));
		}
	}
	selectWindow("Log");
	saveAs("Text", myDir+"LSF_Result_COMBO_RECO.txt");
	selectWindow("Log");
	run("Close");
}

//Printing of the FINAL REPORT
if(context == 1) //MAMMO
{
	print("MAMMO CONTEXT:");
	//AEC
	print("AEC Result");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("\n\nmm\tMODE\t ANODE\tFILTER\tkV\tmAs\tExp. Time (ms)\tAnodic Current (mA)\tESAK (mGy)\tHVL (mmAl)\tMGD (mGy)\tDose Rate (mGy/s)");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1)//MAMMO
		{
			print(d2s(mm_array[j],0),"\t"+Exp_mode+"\t"+Anode+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(exp_time_array[j],0)+"\t"+d2s(mA_array[j],0)+"\t"+d2s(ESAK_array[j],2)+"\t"+d2s(HVL_array[j],3)+"\t"+d2s(MGD_array[j]*100,2)+"\t"+d2s((ESAK_array[j])/((exp_time_array[j])/1000),2));
		}
	}
	if (result_exp_time == "Pass")
	{
		print("\nTube Efficiency OK. Exposure time for 40 mm is under the limit of 1s");
	}
	if (result_exp_time == "Fail")
	{
		print("\nTube Efficiency Fail. Exposure time for 40 mm is over the limit of 1s");
	}
	//SNR-CNR
	print("\n****************************************************************************************");
	print("\nSNR-CNR Result");
	print("\nmm\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tLSF\tSNR Result");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1)//MAMMO
		{
		   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(MeanBKG_array[j],1)+"\t"+d2s(StdBKG_array[j],2)+"\t"+d2s(SNR_BKG_array[j],1)+"\t"+d2s(MeanDetail_array[j],1)+"\t"+d2s(StdDetail_array[j],2)+"\t"+d2s(SNR_Object_array[j],1)+"\t"+d2s(Contrast_array[j],3)+"\t"+d2s(CNR_array[j],1)+"\t"+d2s(lsf_array[j],3)+"\t"+result_SNR_array[j]);
		}
	}
	print("\nTest Limit SNR > 60");
	
	//Homogeneity
	print("\n****************************************************************************************");
	print("\nHOMOGENEITY Result");	
	
	print("\nmm\tFILTER\tkV\tmAs\tROI1\tROI2\tROI3\tROI4\tROI5\tAverage\tDev max % Mean");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1) //MAMMO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(MeanROI1[j],1)+"\t"+d2s(MeanROI2[j],2)+"\t"+d2s(MeanROI3[j],1)+"\t"+d2s(MeanROI4[j],1)+"\t"+d2s(MeanROI5[j],2)+"\t"+d2s(MeanTot[j],1)+"\t"+d2s(Max_Mean_Dev[j],3));
		}
	}
	print("\nmm\tFILTER\tkV\tmAs\tSNR1\tSNR2\tSNR3\tSNR4\tSNR5\tAverage SNR\tDev % SNR");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1) //MAMMO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(SNR_ROI1[j],1)+"\t"+d2s(SNR_ROI2[j],1)+"\t"+d2s(SNR_ROI3[j],1)+"\t"+d2s(SNR_ROI4[j],1)+"\t"+d2s(SNR_ROI5[j],1)+"\t"+d2s(SNRTot[j],1)+"\t"+d2s(Max_SNR_Dev[j],3));
		}
	}
	print("\nHomogeneity Test Result for Mean Grey Level\t"+result_Homogeneity_Mean);
	print("Limit Dev % max for average grey level < 5%");
	print("\nHomogeneity Test Result for SNR\t"+result_Homogeneity_SNR);
	print("Limit Dev % max for SNR < 15%");
	
	//Grid Visibility
	print("\n****************************************************************************************");
	print("\nGrid Visibility Result");	
	print("\nmm\tFILTER\tkV\tmAs\tExp. Time (ms)\tGrid Peak\tNPS Ground\tPeak-Ground Ratio");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1)//MAMMO
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(exp_time_array[j],0)+"\t"+d2s(Grid_Peak_array[j],-1)+"\t"+d2s(NPS_Ground_array[j],-1)+"\t"+d2s(grid_visibility_array[j],1));
		}
	}
	print("\nLimit Peak-Ground Ratio < 20");
	
	//LSF-MTF
	print("\n****************************************************************************************");
	print("\nLSF-MTF Result");
	print("\nmm\tFILTER\tkV\tmAs\tangle (degree)\tFHWM (mm)");
	for(j = 0; j < list.length; j++)
	{
	   print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(wire_angle_array[j],3)+"\t"+d2s(lsf_array[j],3));
	}
	print("\nTest Limit FWHM < 0.150 mm");
	selectWindow("Log");
	saveAs("Text", myDir+"Summary of Result - MAMMO MODE.txt");
	selectWindow("Log");
	run("Close");
}

if(context == 2) //TOMO
{
	print("TOMO CONTEXT:");
	//AEC-CNR TOMO
	print("\nAEC Result:");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	
	print("\nAEC-CNR PROJECTION:");
	print("\nmm\tAngle (°)\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tSNR Result\t");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 2)//PROJECTION
		{
			print(d2s(mm_array[j],0)+"\t"+d2s(angle_array[j],2)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(uAs_array[j]/1000,2)+"\t"+d2s(MeanBKG_array[j],1)+"\t"+d2s(StdBKG_array[j],2)+"\t"+d2s(SNR_BKG_array[j],1)+"\t"+d2s(MeanDetail_array[j],1)+"\t"+d2s(StdDetail_array[j],2)+"\t"+d2s(SNR_Object_array[j],1)+"\t"+d2s(Contrast_array[j],3)+"\t"+d2s(CNR_array[j],2)+"\t"+result_SNR_array[j]);
		}
	}
	print("\nTest Limit SNR > 25");
	
	print("\nAEC-CNR RECONSTRUCTIONS:");
	print("\nmm\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tsigma LSF\tSNR Result\tCNR ASF");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3)//RECO
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(BackGround_mean_array[j],1)+"\t"+d2s(BackGround_std_array[j],2)+"\t"+d2s(BackGround_SNR_array[j],1)+"\t"+d2s(Al_mean_array[j],1)+"\t"+d2s(Al_std_array[j],2)+"\t"+d2s(Al_SNR_array[j],1)+"\t"+d2s(Contrast_rec_array[j],3)+"\t"+d2s(CNR_asf_focused_slice[j],2)+"\t"+d2s(lsf_array[j],3)+"\t"+result_SNR_array[j]+"\t"+d2s(CNR_rec_array[j],2));
		}
	}
	print("\nTest Limit SNR > 50");
	
	//HOMOGENEITY TOMO
	print("\n****************************************************************************************");
	print("\nHOMOGENEITY RECONSTRUCTIONS:");
	print("\nmm PMMA\t FILTER\t kV\t mAs tot.\t ROI1\t ROI2\t ROI3\t ROI4\t ROI5\t Average\t Dev max % Mean");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3) //RECO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(MeanROI1[j],1)+"\t"+d2s(MeanROI2[j],2)+"\t"+d2s(MeanROI3[j],1)+"\t"+d2s(MeanROI4[j],1)+"\t"+d2s(MeanROI5[j],2)+"\t"+d2s(MeanTot[j],1)+"\t"+d2s(Max_Mean_Dev[j],3));
		}
	}
	print("\nmm\tFILTER\tkV\tmAs\tSNR1\tSNR2\tSNR3\tSNR4\tSNR5\tAverage SNR\tDev % SNR");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3) //RECO Average SNR Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(SNR_ROI1[j],1)+"\t"+d2s(SNR_ROI2[j],1)+"\t"+d2s(SNR_ROI3[j],1)+"\t"+d2s(SNR_ROI4[j],1)+"\t"+d2s(SNR_ROI5[j],1)+"\t"+d2s(SNRTot[j],1)+"\t"+d2s(Max_SNR_Dev[j],3));
		}
	}
	print("\nHomogeneity Test Result for Mean Grey Level\t"+result_Homogeneity_Mean);
	print("Limit Dev % max for average grey level < 5%");
	print("\nHomogeneity Test Result for SNR\t"+result_Homogeneity_SNR);
	print("Limit Dev % max for SNR < 15%");
	
	//LSF-MTF TOMO
	print("\n****************************************************************************************");
	print("\nLSF-MTF RECONSTRUCTIONS:");
	if(LSF_cal == 1) //TOMO
	{
		print("\nmmPMMA\tFILTER\tkV\tmAstot.\tangle (degree)\tFHWM (mm)");
		for(j = 0; j < list.length; j++)
		{
			if (ImageType_array[j] == 3)//RECO
			{
				print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(wire_angle_array[j],3)+"\t"+d2s(lsf_array[j],3));
			}
		}
	}
	print("\nTest Limit FWHM < 0.250 mm");
	selectWindow("Log");
	saveAs("Text", myDir+"Summary of Result - TOMO MODE.txt");
	selectWindow("Log");
	run("Close");
}


if(context == 3) //COMBO
{
	print("COMBO CONTEXT:");
	//AEC-CNR COMBO
	print("\nAEC Result:");
	print("Images size:\t"+d2s(nWidth,0)+"\t x\t"+d2s(nHeight,0)+"\t"+d2s(BitStored,0)+"Bit");
	print("Machine Model:\t"+MachineModel);
	print("Machine Serial Number:\t"+d2s(GiottoSerialNumber,0));
	print("Detector ID:\t"+DetectorID);
	print("SoftwareVersion:\t"+SoftwareVersion);
	print("Last Calibration Date:\t"+d2s(LastCalibDate,0));
	print("SID:\t"+d2s(SID,1));
	print("SOD:\t"+d2s(SOD,1));
	print("SAD:\t"+d2s(SAD,1));
	print("FSP_X:\t"+d2s(FSP_X,0));
	print("FSP_Y:\t"+d2s(FSP_Y,0));
	
	print("\nAEC-CNR COMBO MAMMO:");
	print("\nmm\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tLSF\tSNR Result");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1)//MAMMO
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(MeanBKG_array[j],1)+"\t"+d2s(StdBKG_array[j],2)+"\t"+d2s(SNR_BKG_array[j],1)+"\t"+d2s(MeanDetail_array[j],1)+"\t"+d2s(StdDetail_array[j],2)+"\t"+d2s(SNR_Object_array[j],1)+"\t"+d2s(Contrast_array[j],3)+"\t"+d2s(CNR_array[j],1)+"\t"+d2s(lsf_array[j],3)+"\t"+result_SNR_array[j]);
		}
	}
	print("\nTest Limit SNR > 25");
	
	print("\nAEC-CNR PROJECTION:");
	print("\nmm\tAngle (°)\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tSNR Result\t");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 2)//PROJECTION
		{
			print(d2s(mm_array[j],0)+"\t"+d2s(angle_array[j],2)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(uAs_array[j]/1000,2)+"\t"+d2s(MeanBKG_array[j],1)+"\t"+d2s(StdBKG_array[j],2)+"\t"+d2s(SNR_BKG_array[j],1)+"\t"+d2s(MeanDetail_array[j],1)+"\t"+d2s(StdDetail_array[j],2)+"\t"+d2s(SNR_Object_array[j],1)+"\t"+d2s(Contrast_array[j],3)+"\t"+d2s(CNR_array[j],2)+"\t"+result_SNR_array[j]);
		}
	}
	print("\nTest Limit SNR > 25");
	
	print("\nAEC-CNR RECONSTRUCTIONS:");
	print("\nmm\tFILTER\tkV\tmAs\tP.V. PMMA\ts PMMA\tSNR PMMA\tP.V. Al\ts Al\tSNR Al\tContrast\tCNR\tsigma LSF\tSNR Result\tCNR ASF");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3)//RECO
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(BackGround_mean_array[j],1)+"\t"+d2s(BackGround_std_array[j],2)+"\t"+d2s(BackGround_SNR_array[j],1)+"\t"+d2s(Al_mean_array[j],1)+"\t"+d2s(Al_std_array[j],2)+"\t"+d2s(Al_SNR_array[j],1)+"\t"+d2s(Contrast_rec_array[j],3)+"\t"+d2s(CNR_asf_focused_slice[j],2)+"\t"+d2s(lsf_array[j],3)+"\t"+result_SNR_array[j]+"\t"+d2s(CNR_rec_array[j],2));
		}
	}
	print("\nTest Limit SNR > 50");
	
	//HOMOGENEITY COMBO
	print("\n****************************************************************************************");
	print("\nHOMOGENEITY COMBO MAMMO:");
	print("\nmm PMMA\t FILTER\t kV\t mAs\t ROI1\t ROI2\t ROI3\t ROI4\t ROI5\t Average\t Dev max % Mean");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1) //MAMMO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(MeanROI1[j],1)+"\t"+d2s(MeanROI2[j],2)+"\t"+d2s(MeanROI3[j],1)+"\t"+d2s(MeanROI4[j],1)+"\t"+d2s(MeanROI5[j],2)+"\t"+d2s(MeanTot[j],1)+"\t"+d2s(Max_Mean_Dev[j],3));
		}
	}
	print("\nmm\tFILTER\tkV\tmAs\tSNR1\tSNR2\tSNR3\tSNR4\tSNR5\tAverage SNR\tDev % SNR");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 1) //MAMMO Average SNR Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(SNR_ROI1[j],1)+"\t"+d2s(SNR_ROI2[j],1)+"\t"+d2s(SNR_ROI3[j],1)+"\t"+d2s(SNR_ROI4[j],1)+"\t"+d2s(SNR_ROI5[j],1)+"\t"+d2s(SNRTot[j],1)+"\t"+d2s(Max_SNR_Dev[j],3));
		}
	}
	print("\n****************************************************************************************");
	print("\nHOMOGENEITY RECONSTRUCTIONS:");
	print("\nmm\tFILTER\tkV\tmAs tot.\tROI1\tROI2\tROI3\tROI4\tROI5\tAverage\tDev max % Mean");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3) //RECO Average Pixel Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(MeanROI1[j],1)+"\t"+d2s(MeanROI2[j],2)+"\t"+d2s(MeanROI3[j],1)+"\t"+d2s(MeanROI4[j],1)+"\t"+d2s(MeanROI5[j],2)+"\t"+d2s(MeanTot[j],1)+"\t"+d2s(Max_Mean_Dev[j],3));
		}
	}
	print("\n mm\tFILTER\tkV\tmAs\tSNR1\tSNR2\tSNR3\tSNR4\tSNR5\tAverage SNR\tDev % SNR");
	for(j = 0; j < list.length; j++)
	{
		if (ImageType_array[j] == 3) //RECO Average SNR Value
		{
			print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(SNR_ROI1[j],1)+"\t"+d2s(SNR_ROI2[j],1)+"\t"+d2s(SNR_ROI3[j],1)+"\t"+d2s(SNR_ROI4[j],1)+"\t"+d2s(SNR_ROI5[j],1)+"\t"+d2s(SNRTot[j],1)+"\t"+d2s(Max_SNR_Dev[j],3));
		}
	}
	print("\nHomogeneity Test Result for Mean Grey Level\t"+result_Homogeneity_Mean);
	print("Limit Dev % max for average grey level < 5%");
	print("\nHomogeneity Test Result for SNR\t"+result_Homogeneity_SNR);
	print("Limit Dev % max for SNR < 15%");
	
	//LSF-MTF COMBO
	print("\n****************************************************************************************");
	print("\nLSF-MTF MAMMO COMBO:");
	if(LSF_cal == 1)
	{
		print("\nmm\tFILTER\tkV\tmAstot.\tangle (degree)\tFHWM (mm)");
		for(j = 0; j < list.length; j++)
		{
			if (ImageType_array[j] == 1)//MAMMO
			{
				print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_array[j],0)+"\t"+d2s(wire_angle_array[j],3)+"\t"+d2s(lsf_array[j],3));
			}
		}
	}
	print("\nTest Limit FWHM < 0.150 mm");
	print("\n****************************************************************************************");
	print("\nLSF-MTF RECONSTRUCTIONS:");
	if(LSF_cal == 1)
	{
		print("\nmm\tFILTER\tkV\tmAstot.\tangle (degree)\tFHWM (mm)");
		for(j = 0; j < list.length; j++)
		{
			if (ImageType_array[j] == 3)//RECO
			{
				print(d2s(mm_array[j],0)+"\t"+Filter_array[j]+"\t"+d2s(kV_array[j],0)+"\t"+d2s(mAs_tot_array[j],0)+"\t"+d2s(wire_angle_array[j],3)+"\t"+d2s(lsf_array[j],3));
			}
		}
	}
	print("\nTest Limit FWHM < 0.250 mm");
	selectWindow("Log");
	saveAs("Text", myDir+"Summary of Result - COMBO MODE.txt");
	selectWindow("Log");
	run("Close");
}
//FUNCTIONS
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
//FFT

//Fourier Transform
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
		theta = isign*(6.28318530717959/mmax); 		//Initialize the trigonometric recurrence.
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
