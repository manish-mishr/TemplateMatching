#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <DrawText.h>

using namespace std;

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedSymbol class may be helpful!
//  Feel free to modify.
//
typedef enum {NOTEHEAD=0, QUARTERREST=1, EIGHTHREST=2} Type;
class DetectedSymbol {
public:
  int row, col, width, height;
  Type type;
  char pitch;
  double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<struct DetectedSymbol> &symbols)
{
  ofstream ofs(filename.c_str());

  for(int i=0; i<symbols.size(); i++)
    {
      const DetectedSymbol &s = symbols[i];
      ofs << s.row << " " << s.col << " " << s.width << " " << s.height << " ";
      if(s.type == NOTEHEAD)
	ofs << "filled_note " << s.pitch;
      else if(s.type == EIGHTHREST)
	ofs << "eighth_rest _";
      else 
	ofs << "quarter_rest _";
      ofs << " " << s.confidence << endl;
    }
}

// Function that outputs a visualization of detected symbols
void  write_detection_image(const string &filename, const vector<DetectedSymbol> &symbols, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];
  for(int i=0; i<3; i++)
    output_planes[i] = input;

  for(int i=0; i<symbols.size(); i++)
    {
      const DetectedSymbol &s = symbols[i];

      overlay_rectangle(output_planes[s.type], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 255, 2);
      overlay_rectangle(output_planes[(s.type+1) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);
      overlay_rectangle(output_planes[(s.type+2) % 3], s.row, s.col, s.row+s.height-1, s.col+s.width-1, 0, 2);

      if(s.type == NOTEHEAD)
	{
	  char str[] = {s.pitch, 0};
	  draw_text(output_planes[0], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[1], str, s.row, s.col+s.width+1, 0, 2);
	  draw_text(output_planes[2], str, s.row, s.col+s.width+1, 0, 2);
	}
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

void print(SDoublePlane img)
{
	  for(int i=0;i<img.rows();i++)
	  {
		  for(int j=0;j<img.cols();j++)
		  {
			  cout<<img[i][j]<<" ";
		  }
		  cout<<"\n";
	  }
	  cout<<"\n\n";
}

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  SDoublePlane output1(input.rows(), input.cols());

  // Convolution code here
  for(int row = 0; row< input.rows(); ++row)
  {
	  for(int col=0; col<input.cols();++col)
	  {
		  double sum=0;
		  for(int k = row_filter.cols()-1, u=col-row_filter.cols()/2; k>=0, u<=col+row_filter.cols()/2; --k, ++u)
		  {
			  sum += u>=0 && u<input.cols() ? input[row][u] * row_filter[0][k] : input[row][col]* row_filter[0][k];
		  }

		  output1[row][col] = sum;
	  }
  }

  SDoublePlane output2(input.rows(), input.cols());
  for(int col=0; col<input.cols(); ++col)
  {
	  for(int row=0; row<input.rows(); ++row)
	  {
		  double sum=0;
		  for(int l = col_filter.rows()-1, v = row - col_filter.rows()/2; l>=0, v<=row+col_filter.rows()/2; --l, ++v)
		  {
			  sum += v>=0 && v<output1.rows()? output1[v][col]*col_filter[l][0]: output1[row][col]*col_filter[l][0];
		  }

		  output2[row][col] = sum;
	  }
  }
  return output2;
}

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
  SDoublePlane output(input.rows(), input.cols());

  for(int i=0; i<input.rows();++i)
  {
	  for(int j=0;j<input.cols();++j)
	  {
		  double sum=0;
		  for(int k=filter.rows()-1, u=i-filter.rows()/2; k>=0, u<= i+filter.rows()/2; --k, ++u)
		  {
			  for(int l=filter.cols()-1, v=j-filter.cols()/2; l>=0, v<= j+filter.cols()/2; --l, ++v)
			  {
				  int p = u, q = v;
				  if(p<0 || p >=input.rows()) p = i;
				  if(q<0 || q>input.cols()) q = j;

				  sum+=input[p][q]*filter[k][l];
			  }
		  }
		  output[i][j]= sum;
	  }
  }

  // Convolution code here
  
  return output;
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement a sobel gradient estimation filter with 1-d filters
  

  return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
  SDoublePlane output(input.rows(), input.cols());

  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
  
  return output;
}

void test_colvolution()
{
	SDoublePlane img(5,5);
	for(int i=1; i<4;i++)
		for (int j=1;j<4;j++)
			img[i][j]=16;
	print(img);

	SDoublePlane row_filter(1,3), col_filter(3,1),filter(3,3);
	row_filter[0][0]=.25; row_filter[0][1]=.5; row_filter[0][2]=.25;
	col_filter[0][0]=.25; col_filter[1][0]=.5; col_filter[2][0]=.25;
	filter[0][0]=1/16.0; filter[0][1]=1/8.0; filter[0][2]=1/16.0;
	filter[1][0]=1/8.0; filter[1][1]=1/4.0; filter[1][2]=1/8.0;
	filter[2][0]=1/16.0; filter[2][1]=1/8.0; filter[2][2]=1/16.0;

	SDoublePlane output = convolve_separable(img,row_filter,col_filter);
	print(output);

	SDoublePlane output2 = convolve_general(img,filter);
	print(output2);

}

//
// This main file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//
int main(int argc, char *argv[])
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
  
  // test step 2 by applying mean filters to the input image
  SDoublePlane mean_filter(3,3);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      mean_filter[i][j] = 1/9.0;
  SDoublePlane output_image = convolve_general(input_image, mean_filter);
  SImageIO::write_png_file("mean_filtered.png",output_image,output_image,output_image);

  test_colvolution();


  // randomly generate some detected symbols -- you'll want to replace this
  //  with your symbol detection code obviously!
  vector<DetectedSymbol> symbols;
  for(int i=0; i<10; i++)
    {
      DetectedSymbol s;
      s.row = rand() % input_image.rows();
      s.col = rand() % input_image.cols();
      s.width = 20;
      s.height = 20;
      s.type = (Type) (rand() % 3);
      s.confidence = rand();
      s.pitch = (rand() % 7) + 'A';
      symbols.push_back(s);
    }

  write_detection_txt("detected.txt", symbols);
  write_detection_image("detected.png", symbols, input_image);
}
