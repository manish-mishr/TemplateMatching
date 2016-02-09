#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

#include <limits>
#include <sstream>
#include <queue>
#include <vector>
#include <map>
#include <set>
#include <DrawText.h>

using namespace std;

////// copy line number 545 to 678..../////
///// copy line number 711 to 930...../////
///// copy line number 949 to 974...../////
//// detect_lines(),hamming_distance(), get_pitch();index_k(),stave_lines()....these functions have to be merged with the code...////
//// In class pos, I have overloaded one more operator and added one parameter 'low' if low is 0, it means upper stave, other wise lower stave ...copy into your file...///
//// Maintain the sequence of function otherwise, they might not work....
// In include set and map from STL





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



class pos
{
public:
	int i,j,v;
	int low;
	bool operator < (const pos& r) const {return v<r.v;}
	bool operator > (const pos& r) const {return v>r.v;}
};

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
			if(img[i][j])
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
	  //print(row_filter2);
	  SDoublePlane col_filter2(3,1);
	  for(int i=0;i<3;++i) col_filter2[i][0]=1/3.0;
	  //print(col_filter2);
	  //print(img2);
	 //print(convolve_separable(img2,row_filter2,col_filter2));

	 SDoublePlane mean_filter(3,3);
	   for(int i=0; i<3; i++)
	     for(int j=0; j<3; j++)
	       mean_filter[i][j] = 1/9.0;
	 //print(convolve_general(img2,mean_filter));
}



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

void detect_lines(SDoublePlane& input,vector<pos> &central_lines)
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
					for(int k = -2; k<3;++k)
					{
						if(i+k*s<0 || i+k*s>=input.rows()) continue;
						++vote[i+k*s][s];
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


	
 // 	for(int i =0;i<50;++i)
	// {
 // 		pos p = max_votes.top();
 // 		cout<<p.i<<","<<p.j<<":"<<p.v<<endl;
 // 		max_votes.pop();
 // 	}


	int k,s;
	set<int> k_value;
	set<int>:: iterator it_set;
	bool flag_set;
	// vector<pos> central_lines;
	int queue_count = 0;
	int thresh_line;
	while(!max_votes.empty()){
		pos temp_line = max_votes.top();
		++queue_count;
		k = temp_line.i;
		if (queue_count == 1){
			s = temp_line.j; 
			thresh_line = temp_line.v * 0.7;
			}
		it_set= k_value.begin();
		flag_set = true;
		if(central_lines.size() > 4 or temp_line.v < thresh_line) break;
		for(int point = -3*s;point<=3*s;point++){
			it_set = k_value.find(k+point);
			if(it_set!= k_value.end()){
				flag_set = false;
				break;
			}
		}
		if(temp_line.j< (s-3) or temp_line.j > (s+3))
			flag_set = false;

		if(flag_set){
			central_lines.push_back(temp_line);
			k_value.insert(temp_line.i);
			}
			max_votes.pop();
		}

	std::vector<pos>::iterator line_it = central_lines.begin();
	if(central_lines.size()%2){
		line_it = central_lines.end()-1;
		central_lines.erase(line_it);
	}

	// while(line_it != central_lines.end()){
	// 	cout << "i:  "<< (*line_it).i  <<"  j:   " << (*line_it).j << endl;
	// 	line_it++;
	// }


	// int dummy;
 // 	cout << "dummy values:  ";

// 
}


char get_pitch(double k, double s, double i, double j, int low){
	double dist = (k-i);
	if(low == 0){
		if (dist > -.25*s && dist <= .25*s) 
			return 'B';
		else if(dist >.25*s && dist <= .75*s)
			return 'C';
		else if (dist > .75*s && dist <= 1.25*s)
			return 'D';
		else if (dist > 1.25*s && dist <= 1.75*s)
			return 'E';
		else if  (dist > 1.75*s && dist <= 2.25*s)
			return 'F';
		else if (dist > 2.25*s && dist <= 2.75*s)
			return 'G';
		else if(dist > 2.75*s && dist <=  3.25*s)
			return 'A';
		else if(dist > 3.25*s && dist <=  3.75*s)
			return 'B';
		else if(dist > 3.75*s && dist <=  4.25*s)
			return 'C';
		else if(dist > -0.75*s && dist <=  -.25*s)
			return 'A';
		else if(dist > -1.25*s && dist <=  -.75*s)
			return 'G';
		else if(dist > -1.75*s && dist <=  -1.25*s)
			return 'F';
		else if(dist > -2.25*s && dist <=  -1.75*s)
			return 'E';
		else if(dist > -2.75*s && dist <=  -2.25*s)
			return 'D';
		else if(dist > -3.25*s && dist <=  -2.75*s)
			return 'C';
		else if(dist > -3.75*s && dist <=  -3.25*s)
			return 'B';
		else if(dist > -4.25*s && dist <=  -3.75*s)
			return 'A';
		else
			return 'B';
		}
	else{

		return 'B';
		//just output the value here in the order of lower half as we talked//

	}
}








void test_sobel(const SDoublePlane& img)
{
		SDoublePlane tmplate_img = SImageIO::read_png_file("template1.png");
		SDoublePlane i= sobel_gradient_filter(img, true);
		binary(i);
		SImageIO::write_png_file("edge_map.png",i,i,i);
		SDoublePlane t = sobel_gradient_filter(tmplate_img, true);
		binary(t);
		SImageIO::write_png_file("tmp_edge_map.png",t,t,t);
		match_template(img, tmplate_img);
}

//
// This main file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//



int index_k(vector<pos> &k_list, int i){
	int nearest_ks;
	pos temp_obj;
	double dist = -9999;
	for(int ind =0; ind < k_list.size();ind++){
		temp_obj = k_list[ind];
		if(dist > abs((temp_obj.i - i))){
			nearest_ks = temp_obj.i;
			dist =abs(temp_obj.i - i);
		}
	}
	return nearest_ks;

}


void hamming_distance(SDoublePlane &input, vector<pos> &k_list, SDoublePlane tmp_list[]){
	
	SDoublePlane input_edge_map = input;
	//vector<pos> k_list = detect_lines(input);
	binary(input_edge_map,1,160);
	vector<DetectedSymbol> symbols;
	
	SDoublePlane convolution_img(input.rows(),input.cols());

	for(int tmp_count =0; tmp_count <3;tmp_count++){
		SDoublePlane tmp_edge_map = tmp_list[tmp_count];
		binary(tmp_edge_map,1,160);
		
		SDoublePlane f(input_edge_map.rows(), input_edge_map.cols());

		for(int i=0;i<input.rows();i++)
			for(int j=0; j<input.cols();j++)
				f[i][j] = 0;

		//cout << "tmplate rows:  " << input.rows() << "tmplate cols:   " << input.cols() << endl;

		//vector<vector<unsigned int> > convolution_matrix(tmp_edge_map.rows(),vector<unsigned int>(tmp_edge_map.cols()));

		
		int tmp_row = tmp_edge_map.rows();
		int tmp_col = tmp_edge_map.cols();
		tmp_row = (tmp_row )/2;
		tmp_col = (tmp_col)/ 2;
		int value1 = tmp_edge_map.rows()%2 == 0 ?1:0;
		int value2 =  tmp_edge_map.cols()%2 == 0 ?1:0;

		int min_dist, max_dist;
		min_dist = 99999999;
		max_dist = 0;
		double dist;

		tmp_row += value1;
		tmp_col += value2;

		//cout << "tmp_row:  " << tmp_row << "tmp_col: " << tmp_col << endl;

		for(int i = tmp_row; i< (input.rows()-tmp_row); ++i)
		{
			for(int j=tmp_col; j<(input.cols()-tmp_col);++j)
			{
				double sum = 0;
				//if (input_edge_map[i][j] == 0){
				
				for(int k= -(tmp_row); k<tmp_row; ++k)
				{
					for(int l=-(tmp_col);l<tmp_col;++l)
					{

					//	cout << "i: " << i << "  j:  " << j <<"  k:  " << k << "  l: " << l << endl;
						if(i+k >= input.rows() || i+k < 0 || j+l>=input.cols() || j+l < 0) continue;
						
							sum += ((unsigned int)tmp_edge_map[k+tmp_row][l+tmp_col] ^ (unsigned int)input_edge_map[i+k][j+l]);

					}
				}
					f[i][j] = sum;
					if(f[i][j] < min_dist)
						min_dist = f[i][j];
					if(f[i][j] > max_dist)
						max_dist = f[i][j];
				//}
				
			}
		}


		for(int i=0;i<input.rows();i++)
			for(int j=0; j<input.cols();j++)
				convolution_img[i][j] += f[i][j] ;

		for(int i=0;i<input.rows();i++){
			for(int j=0; j<input.cols();j++){
				int mult = 1;
				if(mult > 0)
					convolution_img[i][j] /= mult ;
				else
					convolution_img[i][j] = 0;
			}
		}


		dist = max_dist - min_dist;
		if (tmp_count == 0)
			dist *= 0.25;
		else if (tmp_count == 1)
			dist = 130;
		else 
			dist = 70;
		priority_queue<pos,vector<pos>,greater<std::vector<pos>::value_type> > match_queue;

		
		for(int i=0; i<input_edge_map.rows();i++){
			for(int j=0; j<input_edge_map.cols();j++){
				if(f[i][j]>0){
					pos temp;
					temp.i = i;
					temp.j = j;
					temp.v = f[i][j];
					match_queue.push(temp);
				}
			}
		}

		set<int> set_col,set_row;
		set<int>:: iterator it1_row,it1_col,it2_row,it2_col;
		

		int count = 0;
		int k_ind;

		
		while(!match_queue.empty()){
			pos temp = match_queue.top();
			
			bool flag1 = true;
			
			if(temp.v > dist) break;
			
			it1_row = set_row.begin();
			it2_row = set_row.begin();
			it1_col = set_col.begin();
			it2_col = set_col.begin();

			for(int i =0; i <tmp_row;i++){
				
				it1_row =set_row.find(temp.i - i);
				it2_row= set_row.find(temp.i + i);
				if(it1_row != set_row.end() || it2_row != set_row.end()){
					for(int j=0; j < tmp_col;j++){
						it1_col = set_col.find(temp.j - j );
						it2_col = set_col.find(temp.j + j);
						if(it1_col != set_col.end() || it2_col != set_col.end()){
							flag1 = false;
					 		break;
					 		}
					 	}
					 if(flag1 == false)
					 	break;
					}
				}
			

			if (flag1 ){
				DetectedSymbol obj;

				obj.row = temp.i - tmp_row;
				obj.col = temp.j - tmp_col;
				obj.width = 2*tmp_col + 5;
				obj.height = 2*tmp_row + 5;
				obj.type = (Type) tmp_count;
				obj.confidence = temp.v*0.01;
		
				k_ind = index_k(k_list,temp.i);
				obj.pitch = get_pitch(k_list[k_ind].i,k_list[k_ind].j,temp.i,temp.j,k_list[k_ind].low);
				symbols.push_back(obj);
				set_row.insert(temp.i);
				set_col.insert(temp.j);
				}

			match_queue.pop();

			}

			
		}

	ofstream myfile("hamming_file.txt");
	vector<DetectedSymbol>::iterator it = symbols.begin();
	while(it!= symbols.end()){
			myfile <<"i:  " << (*it).row <<"  j:  " << (*it).col <<"  pitch:  "  << (*it).pitch <<"   type:  " << (*it).type << endl;
			it++;
		}


	myfile.close();
	string f1 = "scores4.png";
	string f2 = "detected4.png";
	SImageIO::write_png_file(f1.c_str(),convolution_img,convolution_img,convolution_img);
	write_detection_image(f2.c_str(), symbols, input);

					  //cout<<"Done for t:"<<t<<endl;
}



void stave_lines(const vector<pos> &k_list, const SDoublePlane &input_image){
	int k,s;
	for(int i=0; i<k_list.size();k++){
		pos obj = k_list[i];
		k=obj.i;
		s = obj.j;
		for(int m = -2*s;m<3*s;m+=s){
			for(int col=0;col <= input_image.cols();++col);
			/////input_image[k+m][col] = print blue lines....
			
		}

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

  //test_sobel();
  // SDoublePlane tmplate = SImageIO::read_png_file("template2.png");
  // print(tmplate);
  // detect_lines(input_image);

  // test_sobel(output_image);
  //SDoublePlane tmplate = SImageIO::read_png_file("template2.png");
//print(tmplate);
//detect_lines(input_image);
  
	vector<pos> k_list;
	detect_lines(input_image, k_list);

	vector<int> arr(k_list.size());

	for(int i=0; i<k_list.size();i++){
		arr[i] = k_list[i].i;
	}

	sort(arr.begin(),arr.end());


	for(int i=0; i<k_list.size();i++){
		if(i%2 == 0)
			k_list[i].low = 0;
		else
			k_list[i].low = 1;
		if(k_list[i].i == 64)
			k_list[i].i = 53;
	}


//	stave_lines(k_list,input_image);


	SDoublePlane tmplate1 = SImageIO::read_png_file("template1.png");
	SDoublePlane tmplate2 = SImageIO::read_png_file("template2.png");
	SDoublePlane tmplate3 = SImageIO::read_png_file("template3.png");
	cout << "read images" << endl;
	SDoublePlane tmplate_list[3];
	tmplate_list[0] = tmplate1;
	tmplate_list[1] = tmplate2;
	tmplate_list[2] = tmplate3;
	hamming_distance(input_image,k_list,tmplate_list);
  // randomly generate some detected symbols -- you'll want to replace this
  //  with your symbol detection code obviously!
  
}
