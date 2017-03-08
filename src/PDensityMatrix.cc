////////////////////////////////////////////////////////
//  
// Reads an ASCII-matrix and samples observables
// in 1-, 2- or 3-dimensional space. The sampled
// observables can be written to PParticles via
// batch scripting
//
////////////////////////////////////////////////////////

#include "PDensityMatrix.h"
#include "PChannel.h"


PDensityMatrix::PDensityMatrix() {
    projector      = NULL;
    current_matrix = 0;
    x = makeStaticData()->GetBatchValue("_x");
    y = makeStaticData()->GetBatchValue("_y");
    z = makeStaticData()->GetBatchValue("_z");
}

Bool_t PDensityMatrix::IsBorder(Int_t bin) {
    //checks if there is a grid border between bin-1 and bin
    if ((bin % axes_size[0]) == 0) {
	return kTRUE;
    }
    return kFALSE;
}

Bool_t PDensityMatrix::GetBin(Double_t *x, Int_t *bins) {
    bins[0] = bins[1] = bins[2] = 0;

    for (int i=0; i<dimension; i++) {
	bins[i] = 0;
	for (int j=0;j<axes_size[i];j++) {
	    bins[i] = j;
	    if (x[i] <= (axes[i])[j]) {
		j = axes_size[i];
	    } 
	}
    }
    
    return kTRUE;
}

Bool_t PDensityMatrix::GetBin(Int_t bin, Int_t *bins) {
    bins[0] = bin % axes_size[0];
    if (dimension > 1) 
	bins[1] = ((bin - bins[0]) / axes_size[0]) % axes_size[1];
    if (dimension > 2) 
	bins[2] = ((bin - bins[0] - bins[1]*axes_size[1]) / (axes_size[0] * axes_size[1])) % axes_size[2];
    return kTRUE;
}

Double_t PDensityMatrix::GetBinWidth(Int_t dim, Int_t bin) {
    if (bin <= 0) {
	return ((axes[dim])[1] - (axes[dim])[0]);
    } else if (bin >= axes_size[dim]-1) {
	return ((axes[dim])[axes_size[dim]-1] - (axes[dim])[axes_size[dim]-2]);
    } 
    return ((axes[dim])[bin+1] - (axes[dim])[bin-1]) / 2.0;
}

Int_t PDensityMatrix::GetRandomBin(Int_t num) {
    Int_t lower_bound = -1;
    Int_t upper_bound = matrixsize - 1;    
    Int_t middle_point;

    while ((upper_bound - lower_bound) > 1) {
	middle_point = lower_bound + (upper_bound-lower_bound)/2;

	Double_t frac = 0;

	if (lower_bound < 0)
	    frac = (matrix_integral[num])[middle_point] /
		(matrix_integral[num])[upper_bound];
	else
	    frac = ((matrix_integral[num])[middle_point] - (matrix_integral[num])[lower_bound]) /
		((matrix_integral[num])[upper_bound] - (matrix_integral[num])[lower_bound]);
	
	if (PUtils::sampleFlat() > frac) {
	    lower_bound = middle_point;
	} else {
	    upper_bound = middle_point;
	}
    }
    return lower_bound + 1;
}

Bool_t PDensityMatrix::Do(const char * command) {

    if (!projector) projector = new PProjector();
    return projector->AddCommand(command);

    return kTRUE;
}

Bool_t PDensityMatrix::ReadDensityMatrix(char *filename, Int_t dim, Bool_t use_bin_width,
					 Double_t min_selection, Double_t max_selection) {
    //Reads a density matrix which must be organized in one of the the two following way
    //
    //Method 1:
    //x [y] [z] f1 f2 f3 ....
    //
    //The first column contains the bin center on the x-axis, the following
    //columns the values of 1 or more matrices. Optionally, for 2- and 3-dimensional
    //matrixes, the first colum might be followed by the y- and z-value of the bin center
    //
    //The matrixes support variable bin width, the flag "use_bin_width" should be used
    //to select of the sample statistics should be corrected by the bin width
    //
    //Method 2:
    //"section value 1"
    //x [y] [z] f0 f1 f2 ....
    //...
    //"section value 2"
    //x [y] [z] f0 f1 f2 ....
    //...
    //
    //Here, the syntax is like in method 1, but if lines have only one number, 
    //they are considered to be a "selection number". Which of the sections should be read
    //and used to fill the internal matrix, can be selected by "min_selection" and "max_selection"
    //The section is used if "min_selection" <= "section value" <= "max_selection"
    //
    //In your macro, to select one of the matrixes, use
    //  matrix->SetMatrix(N);
    //where N corresponds to fN, e.g.
    //  matrix->SetMatrix(2);
    //uses the values f2 as defined above

    if (pluto_global::verbosity) {
        Info("ReadDensityMatrix", "Analysing the file...");
    }

    Double_t  numbers[16];
    dimension = dim;

    for (int i=0; i<dim; i++) {
	axes[i]      = new Double_t[DENSITYMATRIX_MAX];
	axes_size[i] = 0;
    }

    Bool_t found_selection      = kFALSE;
    Bool_t found_single_numbers = kFALSE;
    FILE * file                 = fopen(filename,"r");
    Int_t  line_number = 0;

    if (!file) {
	Error("ReadDensityMatrix", "File '%s' not found", filename);
	return kFALSE;
    }
    
    Int_t numargs = 0;
    char line [1024];
    
    while(!found_selection && (fgets (line, sizeof line, file) != NULL)) {
	line_number++;
	numargs = sscanf(line, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
			 &(numbers[0]), &(numbers[1]), &(numbers[2]), &(numbers[3]), &(numbers[4]), 
			 &(numbers[5]), &(numbers[6]), &(numbers[7]), &(numbers[8]), &(numbers[9]), 
			 &(numbers[10]), &(numbers[11]), &(numbers[12]), &(numbers[13]), &(numbers[14]), 
			 &(numbers[15]));
	
	if (numargs == 1) {
	    found_single_numbers = kTRUE;
	    if (numbers[0] >= min_selection && numbers[0] <= max_selection)
		found_selection = kTRUE;
	}
    }

    if (!found_single_numbers) { //no selection at all
	fclose(file);
	file = fopen(filename,"r");
	line_number = 0;
    } 

    Int_t num_matrices = 0;

    //continuing with reading dimensions of matrix
    while(fgets (line, sizeof line, file) != NULL) {
	line_number++;
	numargs = sscanf(line, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
			 &(numbers[0]), &(numbers[1]), &(numbers[2]), &(numbers[3]), &(numbers[4]), 
			 &(numbers[5]), &(numbers[6]), &(numbers[7]), &(numbers[8]), &(numbers[9]), 
			 &(numbers[10]), &(numbers[11]), &(numbers[12]), &(numbers[13]), &(numbers[14]), 
			 &(numbers[15]));
	
	if (numargs == 1) {
	    if (numbers[0] >= min_selection && numbers[0] <= max_selection)
		found_selection = kTRUE;
	    else
		found_selection = kFALSE;
	} else if (found_selection && numargs > 1 ) {
	    if (numargs <= dim) {
		Error("ReadDensityMatrix", "Not enough values found");
		cout << "Line number " << line_number << ": " << line << endl;;
		return kFALSE;
	    }
	    for (int i=0; i<dim; i++) {
		if (!axes_size[i] || numbers[i] > (axes[i])[axes_size[i]-1]) {
		    if (axes_size[i] == DENSITYMATRIX_MAX) {
			Error("ReadDensityMatrix", "Size too small for dimension %i", i);
			return kFALSE;
		    }
		    (axes[i])[axes_size[i]] = numbers[i];
		    axes_size[i]++;		    
		}
	    }
	    if (!num_matrices) {
		num_matrices = numargs - dim;
	    } else if (num_matrices && num_matrices != numargs - dim) {
		Error("ReadDensityMatrix", "Number of matrices do not match, it was %i, and is now %i",
		      num_matrices, numargs - dim);
		cout << "Line number " << line_number << ": " << line << endl;;
		return kFALSE;
	    }
	}
    }

    if (pluto_global::verbosity) {
        if (!found_single_numbers)
            Info("ReadDensityMatrix", "...done (no sections found)");
        else
            Info("ReadDensityMatrix", "...done");

        Info("ReadDensityMatrix", "Dimension 1 (_x) has %i bins within the range [%f,%f]",
             axes_size[0], (axes[0])[0], (axes[0])[axes_size[0]-1]);
        if (dimension > 1)
            Info("ReadDensityMatrix", "Dimension 2 (_y) has %i bins within the range [%f,%f]",
                 axes_size[1], (axes[1])[0], (axes[1])[axes_size[1]-1]);
        if (dimension > 2)
            Info("ReadDensityMatrix", "Dimension 3 (_z) has %i bins within the range [%f,%f]",
                 axes_size[2], (axes[2])[0], (axes[2])[axes_size[2]-1]);
    }

    //TODO: check axis size against dim and DENSITYMATRIX_MAX_MATRICES

    fclose(file);

    Double_t x[3];
    Int_t    bins[3];

    //based on the dimensions, let's create the matrices
    matrixsize = 0;
    for (int i=0; i<num_matrices; i++) {
	if (dim == 1)
	    matrixsize = axes_size[0];
	else if (dim == 2)
	    matrixsize = axes_size[0]*axes_size[1];
	else 
	    matrixsize = axes_size[0]*axes_size[1]*axes_size[2];
	matrix[i] = new Double_t[matrixsize];
	for (int j=0; j<matrixsize; j++) {
	    (matrix[i])[j] = 0;
	}
    }

    //now filling the content
    file = fopen(filename,"r");
    while(fgets (line, sizeof line, file) != NULL) {
        line_number = 0;
        numargs = sscanf(line, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
                         &(numbers[0]), &(numbers[1]), &(numbers[2]), &(numbers[3]), &(numbers[4]),
                &(numbers[5]), &(numbers[6]), &(numbers[7]), &(numbers[8]), &(numbers[9]),
                &(numbers[10]), &(numbers[11]), &(numbers[12]), &(numbers[13]), &(numbers[14]),
                &(numbers[15]));

        if (numargs == 1) {
            if (numbers[0] >= min_selection && numbers[0] <= max_selection) {
                if (pluto_global::verbosity) {
                    Info("ReadDensityMatrix", "Found section [%f], reading...", numbers[0]);
                }
                found_selection = kTRUE;
            } else {
                if (pluto_global::verbosity) {
                    Info("ReadDensityMatrix", "Found section [%f], skipped!", numbers[0]);
                }
                found_selection = kFALSE;
            }
        } else if (found_selection) {
            for (int i=0; i<dim; i++) {
                x[i] = numbers[i];
            }
            GetBin(x, bins);
            for (int i=0; i<num_matrices; i++) {
                (matrix[i])[bins[0] + bins[1]*axes_size[0] + bins[2]*axes_size[1]*axes_size[0]] += numbers[i+dim];
            }
        }
    }

    //next step is to fill the integral
    for (int i=0; i<num_matrices; i++) {
	matrix_integral[i] = new Double_t[matrixsize];
	for (int j=0; j<matrixsize; j++) {
	    Double_t bin_content = (matrix[i])[j];
	    if (use_bin_width) {
		//fold with bin area
		Int_t bins[3];
		GetBin(j, bins);

		if (dimension == 1) 
		    bin_content *= GetBinWidth(0, bins[0]);
		else if (dimension == 2) 
		    bin_content *= GetBinWidth(0, bins[0]) * GetBinWidth(1, bins[1]);
		else
		    bin_content *= GetBinWidth(0, bins[0]) * GetBinWidth(1, bins[1]) * GetBinWidth(2, bins[2]);
	    }
	   
	    if (j)
		(matrix_integral[i])[j] = bin_content + (matrix_integral[i])[j-1];
	    else 
		(matrix_integral[i])[j] = bin_content;
	}	
    }

    return kTRUE;
}

Bool_t PDensityMatrix::Modify(PParticle ** mstack, int *decay_done, int * num, int stacksize) {
    // Read the particles from the defined stack and copy this to the official
    // particle stream

    PEmbeddedParticles::Modify(mstack, decay_done, num, stacksize);
    Int_t bin = GetRandomBin(current_matrix);

    Int_t xb[3];
    if (!GetBin(bin, xb)) return kFALSE;

    Double_t xf[3];
    for (int i=0; i<dimension; i++) {
	xf[i] = (axes[i])[xb[i]];
	xf[i] += GetBinWidth(i, xb[i])*(PUtils::sampleFlat() - 0.5); //avoid binning effects
    }
    
    *x = xf[0];
    if (dimension>1) *y = xf[1];
    if (dimension>2) *z = xf[2];

    projector->Modify(mstack, decay_done, num, stacksize);
			    
    return kTRUE;
}

ClassImp(PDensityMatrix)
