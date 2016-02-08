#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include<math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include<limits>
#include<sstream>
#include<queue>
#include<vector>
#include<map>
#include<ctime>
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
SDoublePlane convolve_separable_DP(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  SDoublePlane output1(input.rows(), input.cols());

  // Convolution code here
  for(int row = 0; row< input.rows(); ++row)
  {
	  double sum = 0;
	  for(int col=0; col<input.cols();++col)
	  {
		if (col == 0)
		{
			for(int k = row_filter.cols()-1, u=col-row_filter.cols()/2; k>=0, u<=col+row_filter.cols()/2; --k, ++u)
			{
				sum += u>=0 && u<input.cols() ? input[row][u] * row_filter[0][k] : input[row][col]* row_filter[0][k];
			}
		}
		else
		{
			int u =  col-row_filter.cols()/2-1;
			sum -=  row_filter[0][row_filter.cols()-1] * (u >= 0 ? input[row][u] : input[row][0]);
			 u =  col+row_filter.cols()/2;
			 sum += row_filter[0][0] * ( u < input.cols() ? input[row][u] : input[row][input.cols()-1]);
		}
		  output1[row][col] = sum;
	  }
  }

  SDoublePlane output2(input.rows(), input.cols());
  for(int col=0; col<input.cols(); ++col)
  {
	  double sum=0;
	  for(int row=0; row<input.rows(); ++row)
	  {
		if (row == 0)
		{
			for(int l = col_filter.rows()-1, v = row - col_filter.rows()/2; l>=0, v<=row+col_filter.rows()/2; --l, ++v)
			{
				sum += v>=0 && v<output1.rows()? output1[v][col]*col_filter[l][0]: output1[row][col]*col_filter[l][0];
			}
		}
		else
		{
			int v = row - col_filter.rows()/2-1;
			sum -= col_filter[col_filter.rows()-1][0] * ( v >=0 ? output1[v][col] : output1[0][col]);
			 v = row + col_filter.rows()/2;
			sum += col_filter[0][0] * (v < output1.rows() ? output1[v][col] : output1[output1.rows()-1][col]);
		}
		  output2[row][col] = sum;
	  }
  }
  return output2;
}
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
  // Implement a sobel gradient estimation filter with 1-d filters
  SDoublePlane col_filter1(3,1);
  col_filter1[0][0] = 1/8.0; col_filter1[1][0] = 2/8.0; col_filter1[2][0] = 1/8.0;
  SDoublePlane row_filter1(1,3);
  row_filter1[0][0] = 1/8.0; row_filter1[0][1] = 0; row_filter1[0][2] = -1/8.0;

  SDoublePlane col_filter2(3,1);
  SDoublePlane row_filter2(1,3);
  row_filter2[0][0]=1/8.0; row_filter2[0][1]=2/8.0; row_filter2[0][2]=1/8.0;
  col_filter2[0][0]=1/8.0; col_filter2[1][0]=0; col_filter2[2][0]=-1/8.0;

  SDoublePlane sobel_x = convolve_separable(input, row_filter1, col_filter1);
  SImageIO::write_png_file("sobel_x.png",sobel_x,sobel_x,sobel_x);
  SDoublePlane sobel_y = (convolve_separable(input, row_filter2, col_filter2));
  SImageIO::write_png_file("sobel_y.png",sobel_y,sobel_y,sobel_y);

  //now combine gradients
  SDoublePlane output(input.rows(),input.cols());
  for(int i=0;i<input.rows();++i)
  {
	  for(int j=0;j<input.cols();++j)
	  {
		  output[i][j] =(abs(sobel_x[i][j])+abs(sobel_y[i][j])); //sqrt( pow(sobel_x[i][j] , 2) + pow(sobel_y[i][j] , 2) );//
	  }
  }
  //print(output);

  return output;
}

double gamma(double in) { return in == 0 ? std::numeric_limits<double>::max() : 0; }

SDoublePlane calculate_D(const SDoublePlane& img)
{
	vector<pair<int,int> > points;

	
	for(int i = 0; i<img.rows(); ++i)
	{
		for(int j=0; j< img.cols(); ++j)
		{
			if(img[i][j]==1)
			{
				points.push_back(std::make_pair(i,j));
			}
		}
	}

	SDoublePlane D(img.rows(),img.cols());
	for(int i = 0; i<D.rows(); ++i)
	{
		for(int j=0; j< D.cols(); ++j)
		{
			D[i][j]=std::numeric_limits<int>::max();
		}
	}

	for(int p=0; p<points.size(); ++p)
	{
		const pair<int,int>& point = points[p];
		for(int i = 0; i<img.rows(); ++i)
		{
			for(int j=0; j<img.cols(); ++j)
			{				
				double d = sqrt(pow(i-point.first,2)+pow(j-point.second,2));
				if(d<D[i][j])
				{
					D[i][j]=d;
				}
			}
		}
	}
	return D;
}



void binary(SDoublePlane& img, int value = 255, int threshold = 2)
{
	
	for(int i=0;i<img.rows();++i)
	{
		for(int j=0;j<img.cols();++j)
		{
			if(abs(img[i][j])> threshold)
			{
				img[i][j] = value;
				//edge_points.push_back(std::make_pair(i,j))
			}
			else
				img[i][j] = 0;
		}
	}
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
  SDoublePlane output = sobel_gradient_filter(input, true);

  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
  binary(output,1,2);
  return output;
}

SDoublePlane match_template(const SDoublePlane& input, const SDoublePlane& tmplate)
{
	SDoublePlane input_edge_map = find_edges(input);
	SDoublePlane tmp_edge_map = find_edges(tmplate);
	cout<<"Edge map found!\n";
	SDoublePlane D = calculate_D(input_edge_map);
	cout<<"D calculated!\n";
	SDoublePlane f(input.rows(), input.cols());

	for(int i = 0; i< D.rows(); ++i)
	{
		for(int j=0; j<D.cols();++j)
		{
			double sum = 0;
			for(int k=0; k<tmplate.rows(); ++k)
			{
				for(int l=0;l<tmplate.cols();++l)
				{
					if(i+k>= D.rows() || j+l>=D.cols()) continue;
					sum += tmp_edge_map[k][l] * D[i+k][j+l];
				}
			}
			f[i][j] = sum;
		}
	}
	cout<<"template matching done!\n";
	//print(f);

	for(int threshold=3;threshold<6;++threshold)
	{
		vector<DetectedSymbol> symbols;
		for(int i = 0; i< f.rows(); ++i)
		{
			for(int j=0; j<f.cols();++j)
			{
				if(f[i][j]<=threshold)
				{
					  DetectedSymbol s;
					  s.row = i;
					  s.col = j;
					  s.width = 20;
					  s.height = 20;
					  s.type = (Type) (rand() % 3);
					  s.confidence = rand();
					  s.pitch = (rand() % 7) + 'A';
					  symbols.push_back(s);

					  ostringstream convert;
					  convert << threshold;
					  string f1 = "detected_"+convert.str()+".txt";
					  string f2 = "detected_"+convert.str()+".png";
					  write_detection_txt(f1.c_str(), symbols);
					  write_detection_image(f2.c_str(), symbols, input);
				}
			}
		}
		 cout<<"Done for t:"<<threshold<<endl;
	}
	cout<<"Image written\n";
}

void test_colvolution()
{
	SDoublePlane img(5,5);
	for(int i=1; i<4;i++)
		for (int j=1;j<4;j++)
			img[i][j]=16;
	//print(img);

	SDoublePlane row_filter(1,3), col_filter(3,1),filter(3,3);
	row_filter[0][0]=.25; row_filter[0][1]=.5; row_filter[0][2]=.25;
	col_filter[0][0]=.25; col_filter[1][0]=.5; col_filter[2][0]=.25;
	filter[0][0]=1/16.0; filter[0][1]=1/8.0; filter[0][2]=1/16.0;
	filter[1][0]=1/8.0; filter[1][1]=1/4.0; filter[1][2]=1/8.0;
	filter[2][0]=1/16.0; filter[2][1]=1/8.0; filter[2][2]=1/16.0;

	//print(convolve_separable(img,row_filter,col_filter));

	SDoublePlane img2(5,5);
	img2[2][2] = 1;
	  SDoublePlane row_filter2(1,3);
	  for(int i=0;i<3;++i) row_filter2[0][i]=1/3.0;
	  print(row_filter2);
	  SDoublePlane col_filter2(3,1);
	  for(int i=0;i<3;++i) col_filter2[i][0]=1/3.0;
	  print(col_filter2);
	  print(img2);
	 print(convolve_separable(img2,row_filter2,col_filter2));

	 SDoublePlane mean_filter(3,3);
	   for(int i=0; i<3; i++)
	     for(int j=0; j<3; j++)
	       mean_filter[i][j] = 1/9.0;
	 print(convolve_general(img2,mean_filter));
}

class pos
{
public:
	int i,j,v;
	bool operator < (const pos& r) const {return v<r.v;}
};

void detect_lines(const SDoublePlane& input, int N_theta, int N_rad, int k)
{
	SDoublePlane edge_map = find_edges(input);
	int u_c = edge_map.rows()/2;
	int v_c = edge_map.cols()/2;
	double r_max = sqrt(pow(u_c,2)+pow(v_c,2));
	double delta_theta = 3.1416/N_theta;
	double delta_rad = 2*r_max/N_rad;

	SDoublePlane acc(N_theta, N_rad);

	for(int i  =0; i<input.rows(); ++i)
	{
		for(int j = 0; j< input.cols(); ++j)
		{
			if(edge_map[i][j])
			{
				int x = i-u_c;
				int y = j - v_c;

				for(int i_theta=0; i_theta<N_theta;++i_theta)
				{
					double theta= i_theta*delta_theta;
					double r = x*cos(theta)+y*sin(theta);
					int i_rad = N_rad/2 + r/delta_rad;
					++acc[i_theta][i_rad];
				}
			}
		}
	}
	priority_queue<pos> max_votes;
	for(int i  =0; i<acc.rows(); ++i)
	{
		for(int j = 0; j<acc.cols(); ++j)
		{
			pos p;
			p.i=i;
			p.j=j;
			p.v=acc[i][j];
			max_votes.push(p);
		}
	}
	for(int i =0;i<10;++i)
	{
		pos p = max_votes.top();
		cout<<p.i<<","<<p.j<<":"<<p.v<<endl;
		max_votes.pop();
	}

}

void detect_lines(const SDoublePlane& input)
{
	SDoublePlane edge_map = find_edges(input);
	SDoublePlane vote(input.rows(), 50); // letting max value of s = 50
	for(int i  =0; i<input.rows(); ++i)
	{
		for(int j = 0; j< input.cols(); ++j)
		{
			if(edge_map[i][j] == 1)
			{
				for( int s =5; s<50;++s)
				{
					for(int k = 0; k<5;++k)
					{
						if(i+k<0 || i+k>=input.rows()) continue;
						++vote[i+k][s];
					}
				}
			}
		}
	}
	cout<<vote.cols()<<endl;
	priority_queue<pos> max_votes;
	for(int i  =0; i<vote.rows(); ++i)
	{
		for(int j = 0; j<vote.cols(); ++j)
		{
			pos p;
			p.i=i;
			p.j=j;
			p.v=vote[i][j];
			max_votes.push(p);
		}
	}
	for(int i =0;i<10;++i)
	{
		pos p = max_votes.top();
		cout<<p.i<<","<<p.j<<":"<<p.v<<endl;
		max_votes.pop();
	}
}

void test_sobel(const SDoublePlane& img)
{
	  	 SDoublePlane filter(3,3); //Sx
		filter[0][0]=-1; filter[0][1]=0; filter[0][2]=1;
		filter[1][0]=-2; filter[1][1]=0; filter[1][2]=2;
		filter[2][0]=-1; filter[2][1]=0; filter[2][2]=1;

		SDoublePlane sobel = sobel_gradient_filter(img, true);
		SImageIO::write_png_file("sobel_both.png",sobel,sobel,sobel);
		//cout<<"\n\nsobel\n\n";
		SDoublePlane tmplate_img = SImageIO::read_png_file("template1.png");
		match_template(img, tmplate_img);
		SDoublePlane tmplate = sobel_gradient_filter(tmplate_img,true);
		SImageIO::write_png_file("template_sobel.png", tmplate,tmplate,tmplate);
		//print(sobel);

		SDoublePlane row_filter(1,3), col_filter(3,1);
			row_filter[0][0]=.25; row_filter[0][1]=.5; row_filter[0][2]=.25;
			col_filter[0][0]=.25; col_filter[1][0]=.5; col_filter[2][0]=.25;

		//make binary version of sobel
		binary(sobel);
		SImageIO::write_png_file("sobel_binary.png",sobel,sobel,sobel);
		binary(tmplate);
		SImageIO::write_png_file("template_binary.png",tmplate,tmplate,tmplate);
}

//
// This main file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//
void t(int input[2][2])
{


	priority_queue<pos> max_votes;
	for(int i  =0; i<2; ++i)
	{
		for(int j = 0; j<2; ++j)
		{
//			map<int,vector<int> > m;
//			m[vote[i][j]].push_back(i);
//			m[vote[i][j]].push_back(j);
//			max_votes.push(m);
			pos p;
			p.i=i;
			p.j=j;
			p.v=input[i][j];
			max_votes.push(p);
		}
	}
	for(int i =0;i<5;++i){
		cout<<max_votes.top().v<<endl;
		max_votes.pop();
	}
}
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

  SDoublePlane row_filter(1,3);
  for(int i=0;i<3;++i) row_filter[0][i]=1/3.0;
  SDoublePlane col_filter(3,1);
  for(int i=0;i<3;++i) col_filter[i][0]=1/3.0;
  SDoublePlane output_image2 = convolve_separable_DP(input_image, row_filter, col_filter);
  SImageIO::write_png_file("mean_filtered2.png",output_image2,output_image2,output_image2);
  //test_colvolution();
  test_sobel(input_image);
  //SDoublePlane tmplate = SImageIO::read_png_file("template2.png");
//print(tmplate);
//detect_lines(input_image);


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
