
#include "abf_data.h"
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>


// TODO proper destructors for everyone and cleanup after main()

/// Construct gradient field object from an ABF-saved file
ABFdata::ABFdata ( const char *gradFileName ) {

    std::ifstream gradFile;
    int		  n;
	char	  hash;	
	double	  xi;

    std::cout << "Opening file " << gradFileName << " for reading\n";
    gradFile.open( gradFileName );
    if (!gradFile.good()) {
	  std::cerr << "Cannot read from file " << gradFileName << ", aborting\n";
	  exit(1);
    }

	gradFile >> hash;
	if ( hash != '#' ) {
	  std::cerr << "Missing \'#\' sign in gradient file\n"; exit(1);
	}
    gradFile >> Nvars;

    std::cout << "Number of variables: " << Nvars << "\n";

    sizes = new int [Nvars];
    blocksizes = new int [Nvars];
    PBC = new int [Nvars];
    widths = new double [Nvars];
    mins = new double [Nvars];

    scalar_dim = 1;	// total is (n1 * n2 * ... * n_Nvars )

    for (int i=0; i<Nvars; i++) {
	  gradFile >> hash;
	  if ( hash != '#' ) {
		std::cerr << "Missing \'#\' sign in gradient file\n"; exit(1);
	  }
	  gradFile >> mins[i] >> widths[i] >> sizes[i] >> PBC[i];  // format is: xiMin dxi Nbins PBCflag
	  std::cout << "min = " << mins[i] << " width = " << widths[i]
	    << " n = " << sizes[i] << " PBC: " << PBC[i] << "\n";

	  if (sizes[i] == 0) {
	    std::cout << "ERROR: size should not be zero!\n";
	    exit(1);
	  }
	  scalar_dim *= sizes[i];
    }
  
	// block sizes, smallest for the last dimension
	blocksizes[Nvars-1] = 1;
    for (int i=Nvars-2; i>=0; i--) {
	  blocksizes[i] = blocksizes[i+1] * sizes[i+1];
	}

    vec_dim = scalar_dim * Nvars;
    std::cout << "Gradient field has length " << vec_dim << "\n";

    gradients = new double [vec_dim];
    deviation = new double [vec_dim];

    for (unsigned int i=0; i<scalar_dim; i++) {
	  for (unsigned int j=0; j<Nvars; j++) {
		// Read and ignore values of the collective variables
		gradFile >> xi;
	  }
	  for (unsigned int j=0; j<Nvars; j++) {
		// Read and store gradients
		gradFile >> gradients[i * Nvars + j];
	  }
    }
    // Could check for end-of-file string here
    gradFile.close(); 

    // for metadynamics
    bias = new double [scalar_dim];
    histogram = new unsigned int [scalar_dim];
    for (unsigned int i=0; i<scalar_dim; i++) {
	  histogram[i] = 0;
	  bias[i] = 0.0;
    }
}

ABFdata::~ABFdata () {
    delete[] sizes;
    delete[] blocksizes;
    delete[] PBC;
    delete[] widths;
    delete[] mins;
    delete[] gradients;
    delete[] deviation;
    delete[] bias;
    delete[] histogram;
}
unsigned int ABFdata::offset ( const int *pos ) {
    unsigned int index = 0; 

    for (int i=0; i<Nvars; i++) {
	// Check for out-of bounds indices here
	  if ( pos[i] < 0 || pos[i] >= sizes[i] ) {
	    std::cerr << "Out-of-range index: " << pos[i] << " for rank " << i << "\n";
	    exit(1);
	  }
	  index += blocksizes[i] * pos[i];
    }
    // we leave the multiplication below for the caller to do
    // we just give the offset for scalar fields
    // index *= Nvars; // Nb of gradient vectors -> nb of array elts
    return index;
}

void ABFdata::write_histogram ( const char *fileName ) {

    std::ofstream os;
    unsigned int index;
    int		*pos, i;

    os.open (fileName);
    if (!os.good()) {
	  std::cerr << "Cannot write to file " << fileName << ", aborting\n";
	  exit(1);
    }
    pos = new int [Nvars];
    for (i=0; i<Nvars; i++) pos[i] = 0;

    for (index = 0; index < scalar_dim; index++) {
	  // Here we do the Euclidian division iteratively
	  for (i=Nvars-1; i>0; i--) {
		  if (pos[i] == sizes[i]) {
			pos[i] = 0;
			pos[i-1]++;
			os << "\n";
		  }
	  }
	  // Now a stupid check:
	  if (index != offset(pos)) {
		  std::cerr << "Wrong position vector at index " << index << "\n";
		  exit(1);
	  }

	  for (i=0; i<Nvars; i++) {
		  os << mins[i] + widths[i] * (pos[i] + 0.5) << " ";
	  }
	  os << histogram[index] << "\n";
	  pos[Nvars-1]++;   // move on to next position
    }	
    os.close();
	delete [] pos;
}


void ABFdata::write_bias ( const char *fileName ) {
// write the opposite of the bias, with global minimum set to 0

    std::ofstream os;
    unsigned int index;
    int		*pos, i;
    double	maxbias;

    os.open (fileName);
    if (!os.good()) {
	  std::cerr << "Cannot write to file " << fileName << ", aborting\n";
	  exit(1);
    }
    pos = new int [Nvars];
    for (i=0; i<Nvars; i++)
	  pos[i] = 0;

	// Set the minimum value to 0 by subtracting each value from the max
    maxbias = bias[0];
    for (index = 0; index < scalar_dim; index++) {
	  if (bias[index] > maxbias)
		maxbias = bias[index]; 
    }

    for (index = 0; index < scalar_dim; index++) {
	  // Here we do the Euclidian division iteratively
	  for (i=Nvars-1; i>0; i--) {
		  if (pos[i] == sizes[i]) {
			pos[i] = 0;
			pos[i-1]++;
			os << "\n";
		  }
	  }
	  // Now a stupid check:
	  if (index != offset(pos)) {
		  std::cerr << "Wrong position vector at index " << index << "\n";
		  exit(1);
	  }

	  for (i=0; i<Nvars; i++) {
		  os << mins[i] + widths[i] * (pos[i] + 0.5) << " ";
	  }
	  os << maxbias - bias[index] << "\n";
	  pos[Nvars-1]++;   // move on to next position
    }	
    os.close();
	delete [] pos;
}


void ABFdata::write_deviation ( const char *fileName ) {
// write the deviation between computed and initial bias
// i.e. "unaccounted for" input gradients

    std::ofstream os;
    unsigned int index;
    int		*pos, i;
	double	*dev;

    os.open (fileName);
    if (!os.good()) {
	  std::cerr << "Cannot write to file " << fileName << ", aborting\n";
	  exit(1);
    }
    pos = new int [Nvars];
    for (i=0; i<Nvars; i++)
	  pos[i] = 0;

	// start at beginning of array
	dev = deviation;

    for (index = 0; index < scalar_dim; index++) {
	  // Here we do the Euclidian division iteratively
	  for (i=Nvars-1; i>0; i--) {
		  if (pos[i] == sizes[i]) {
			pos[i] = 0;
			pos[i-1]++;
			os << "\n";
		  }
	  }

	  for (i=0; i<Nvars; i++) {
		  os << mins[i] + widths[i] * (pos[i] + 0.5) << " ";
	  }
	  for (i=0; i<Nvars; i++) {
		  os << dev[i] << " ";;
	  }
	  os << "\n";

	  pos[Nvars-1]++;   // move on to next position...
	  dev += Nvars;		// ...also in the array
    }	
    os.close();
	delete [] pos;
}