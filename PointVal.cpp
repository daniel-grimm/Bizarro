/*CSS 487 Computer Vision - Final Project
This class creates a PointVal struct which holds the data for 
an OpenCV struct and a double value associated with it.
Additionally a subset of the comparison operators have been overloaded
to increase the usabiliity of this class in the client code.

@author Daniel Grimm
@author Kevin Ulrich*/

//Imports
#include <opencv2/imgproc/imgproc.hpp>

//namespace
using namespace cv;

/*This creates a data structure consisting of an OpenCV Point object
and a double value.*/
struct PointVal
{

	//constructor
	/*Postconditions: A new PointVal object is created
	with the specified Point object and double value.
	@param p : An OpenCV Point object.
	@param d : A double value.*/
	PointVal(Point p, double d, Size s)
	{
		point = p;
		doubleVal = d;
		size = s;
	}

	/*Preconditions: A non-null PointVal object is passed in
	by reference for comparison with another PointVal object.
	Postconditions: If the value of the right hand side is larger
	than the left hand side, true is returned, else false is returned.
	@param p : A non-null PointVal is passed in for comparison.
	@return bool : If the left hand side is larger return false, else
	return true.*/
	bool operator<(PointVal& p)
	{
		if (doubleVal < p.doubleVal)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/*Preconditions: A non-null PointVal object is passed
	into the method for comparison with another PointVal object.
	Postconditions: If the double value of the right hand side is larger
	than the double value of the left hand side, false is returned, else
	true is returned.
	@param p : A non-null PointVal object.
	@return bool : Returns true if the left hand side double value is greater
	than the right hand side double value.*/
	bool operator>(PointVal& p)
	{
		if (doubleVal > p.doubleVal)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	//Global variables
	Point point;	//point in space
	Size size;		//size of the image
	double doubleVal;	//integer value

};