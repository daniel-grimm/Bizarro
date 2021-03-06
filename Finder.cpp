/*CSS 487 Computer Vision - Final Project
This class reads in a database of images and a database of templates,
and finds the number and type of template in the image.

@author Daniel Grimm
@author Kevin Ulrich*/

//Imports
#include <opencv2/core/core.hpp>				//
#include <opencv2/highgui/highgui.hpp>	//
#include <opencv2/imgproc/imgproc.hpp>	//
#include <vector>										//
#include <sstream>									//
#include <fstream>									//
#include <iostream>									//
#include "PointVal.cpp"

//Namespaces
using namespace cv;
using namespace std;

/*Preconditions: An empty vector<Mat>& is passed into the method
for the purpose of containing all of the template images.
Postconditions: The template images are returned as a vector<Mat>&.
@param vector : An empty vector<Mat>& which will contain an array of templates.
@return vector<Mat>& : The array of templates from the database.*/
vector<Mat>& loadTemplateImages(vector<Mat>& vector)
{
	//For the number of templates in the database
	char character = 'a';

	//Read Image into a Mat object
	for (int i = 0; i < 1; i++)
	{
		//Create the name of the file
		string templateName;
		templateName.push_back(character);
		templateName += ".jpg";

		//Read in the image
		Mat templateImage = imread(templateName, CV_LOAD_IMAGE_COLOR);

		//add the image to the vector
		vector.push_back(templateImage);

		//increment the name of the template
		character = character++;
	}
	return vector;
}

/*Preconditions: An empty vector<Mat>& is passed into the method
for the purpose of containing all of the template images.
Postconditions: The template images are returned as a vector<Mat>&.
@param vector : An empty vector<Mat>& which will contain an array of templates.
@return vector<Mat>& : The array of templates from the database.*/
vector<Mat>& loadTemplateMasks(vector<Mat>& vector)
{
	//For the number of templates (masks) in the database
	char character = 'n';

	//Read Image into a Mat object
	for (int i = 0; i < 1; i++)
	{
		//Create the name of the file
		string templateName;
		templateName.push_back(character);
		templateName += ".jpg";

		//Read in the image
		Mat templateMaskImage = imread(templateName, CV_LOAD_IMAGE_GRAYSCALE);

		//add the image to the vector
		vector.push_back(templateMaskImage);

		//increment the name of the mask
		character = character++;
	}
	return vector;
}

/*Preconditions: An empty vector<Mat>& is passed into the method
for the purpose of having the image database loaded into the data structure.
Postconditions: The image database is returned in a vector<Mat>&.
@param vector : An empty vector<Mat>& which will contain an array of images.
@return vector<Mat>& : The database of images.*/
vector<Mat>& loadInputImages(vector<Mat>& vector)
{
	//For the number of images in the database
	for (int i = 0; i < 2; i++)
	{
		//Read Image into a Mat object
		string imageName = "comic" + to_string(i) + ".jpg";
		Mat image = imread(imageName, CV_LOAD_IMAGE_COLOR);

		//add the image to the vector
		vector.push_back(image);
	}
	return vector;
}

/*Preconditions: A non-null image is passed into the method for scaling up.
Postconditions: The same image that was passed in is returned but has been scaled up.
@param templateImage : The image that is going to be scaled.
@return Mat& : The original image but larger.*/
Mat& changeTemplateScale(Mat& templateImage, double scalingFactor)
{
	//resize the image
	resize(templateImage, templateImage, Size(), scalingFactor, scalingFactor, INTER_LINEAR);

	//return the scaled image
	return templateImage;
}

/*Preconditions: A non-null image is passed into the method for rotation.
Postconditions: The same image is returned by has been rotated slightly.
@param templateImage : The image being rotated.
@return Mat& : The rotated image.*/
Mat& changeTemplateRotation(Mat& templateImage, int degreesOfRotation)
{
	//Rotate the template
	Point2d pictureCenter(templateImage.rows / 2.0, templateImage.cols / 2.0);
	Mat rotationMatrix = getRotationMatrix2D(pictureCenter, -degreesOfRotation, 1.0);
	warpAffine(templateImage, templateImage, rotationMatrix, templateImage.size());

	//return the rotated image
	return templateImage;
}

/*Preconditions:
Postconditions:
*/
void drawGreenBox(Mat& inputImage, const PointVal& pv, const string name)
{
	Point point = pv.point;
	Size size = pv.size;

	rectangle(inputImage, point, Point(point.x + size.width, point.y + size.height), Scalar(0, 255, 127), 2);
	imwrite(name, inputImage);
}

/*Preconditions: Non-null values for an input image and template image are
passed into the method. The dimensions of the input image are assumed to be
of equal or larger size than the template.
Postconditions: The number of times the template appears in the image is returned.
@param inputImage : The image being searched over.
@param templateImage : The template that is attempting to be found in the image.
@param bestMatch : The current best match for the template in the image.
@return PointVal&: The point and value of the best match of the template in the image.*/
PointVal& slideTemplateOverImage(Mat& inputImage, Mat& templateImage, Mat& templateMask, PointVal& bestMatch, int& timesUnchanged, bool& changed)
{

	//Create the result matrix
	int resultRows = inputImage.rows - templateImage.rows + 1;
	int resultCols = inputImage.cols - templateImage.cols + 1;
	Mat result(resultRows, resultCols, CV_8UC1);

	//Match the template to the image
	matchTemplate(inputImage, templateImage, result, CV_TM_CCORR_NORMED, templateMask);

	//Get the minimumLocation of the result matrix
	double minimum, maximum;
	Point minLocation, maxLocation;
	minMaxLoc(result, &minimum, &maximum, &minLocation, &maxLocation, Mat());

	//Retrieve the Point and Value of the minimum location
	PointVal testPoint(maxLocation, maximum, Size(templateImage.cols, templateImage.rows));

	//If a better match is found, update the new point
	if (bestMatch < testPoint)
	{
		cout << "updating best match from " << bestMatch.doubleVal << " to " << testPoint.doubleVal << endl;
		cout << "template size: " << templateImage.size << endl;
		bestMatch = testPoint;
		changed = true;
	}

	return bestMatch;//return the best match of the template to the image
}

Mat findEdges(Mat& image) {

	Mat img = image.clone();

	//create grey image
	cvtColor(img, img, COLOR_BGR2GRAY);

	//blur it up
	GaussianBlur(img, img, Size(5, 5), 2.0, 2.0);

	//Sobel it up
	Sobel(img, img, CV_8U, 1, 0);

	return img;

}


/*Preconditions: Non-null values for the input image and template image are passed
in. The input image is assumed to be of the same dimensions or larger than the 
template image.
Parameter bestMatch is assumed to be initialized to Point(0, 0) and DBL_MAX.
Postconditions: The number of occurrences of the template in multiple scales,
and rotations is returned.
@param inputImage : An image that may or may not contain the template.
@param templateImage : The template that will be searched in the image.
@param bestMatch : A PointVal object storing the initial point and double value of the best match for the template.
@param templatesInImage : A vector<vector<int>> that stores the best match for each template for each image.
@param templateNumber : An int representing which template number this image is on.
@return vector<vector<int>> : Each image and the templates that have been found in it.*/
vector< vector<int> >& templateInImage(Mat& inputImage, Mat& templateImage, Mat& templateMask, PointVal& bestMatch, vector< vector<int> >& templatesInImage, int templateNumber, const string name, PointVal& retPointVal)
{
	Mat image = findEdges(inputImage);
	Mat temp = findEdges(templateImage);

	imwrite("edgeimage.jpg", image);

	//initialize an empty pointval
	PointVal test(Point(0, 0), DBL_MIN, Size(0, 0));

	int timesUnchanged = 0;
	bool changed;
	bool stop = false;

	const int max_iterations = 31;
	const double scalingFactor = 0.95;
	const int degreesOfRotation = 90;

	//For the number of possible template scales
	for (int i = 0; i < max_iterations; i++)
	{

		changed = false;

		//For the number of possible template rotations
		for (int j = 0; j < 360 / degreesOfRotation; j++)
		{
			//Find the best match of the template in the image
			test = slideTemplateOverImage(image, temp, templateMask, bestMatch, timesUnchanged, changed);

			//Rotate the image
			temp = changeTemplateRotation(temp, degreesOfRotation);
			templateMask = changeTemplateRotation(templateMask, degreesOfRotation);

			imwrite("template_" + to_string(templateNumber) + "_" + to_string(i) + "_" + to_string(j) + ".jpg", temp);
		}

		if (changed) 
		{
			timesUnchanged = 0;
		}
		else
		{
			timesUnchanged++;
		}


		if (timesUnchanged == int(max_iterations / ceil(scalingFactor * 3))) 
		{
			break;
		}

		//Scale the image
		temp = changeTemplateScale(temp, scalingFactor);
		templateMask = changeTemplateScale(templateMask, scalingFactor);

	}

	drawGreenBox(inputImage, bestMatch, name);
	cout << "Point: " << bestMatch.point << " Value: " << bestMatch.doubleVal << endl;

	//Add the best match to the vector
	templatesInImage[templatesInImage.size() - 1].push_back(templateNumber);

	//Return the number of times the template is in the image
	return templatesInImage;
}

/*Preconditions: An empty vector<string> is passed into the method which
will contain the name of all the templates
Postconditions: The names of all the template images are loaded into the vector
at each index.
@param templateNames : An empty vector of strings.
@return vector<string>& : The list of template names.*/
vector<string>& loadTemplateNames(vector<string>& templateNames)
{
	//open a file to read in the template names
	string title;
	ifstream file("templateNames.txt");

	//If a valid file, continue
	if (file.is_open())
	{
		//While more to read, keep reading
		while (getline(file, title))
		{
			//Add the template name to the file
			templateNames.push_back(title);
		}

		//close the file
		file.close();
	}
	else
	{
		//Unable to open the file
		cout << "Unable to open file templateNames.txt" << endl;
	}

	//return the array of template names
	return templateNames;
}

/*Preconditions: A data structure describing the name of each template
image is passed in as a vector<vector<int>>& object.
Postconditions: The results of the template matching to the image are printed to a text file.
@param vectorOfMatchedTemplates : For each image, the name of the template image as descriged by
an integer index.
@return void*/
void printResultsToFile(vector< vector<int> >& vectorOfMatchedTemplates, vector< PointVal >& resultPointVals)
{
	vector<string> templateNames;
	templateNames = loadTemplateNames(templateNames);

	ofstream out("results.txt");

	for (int i = 0; i < vectorOfMatchedTemplates.size(); i++)
	{
		out << "Image comic" << to_string(i) << ".jpg" << ":" << endl;
		for (int j = 0; j < vectorOfMatchedTemplates[i].size(); j++)
		{
			out << templateNames[vectorOfMatchedTemplates[i][j]] << " - matched at: " << resultPointVals[j].point << " with certainty: " << resultPointVals[j].doubleVal << endl;
		}
		out << endl;
	}
}

/*Preconditions: This program must be run from the command line.
Postconditions: The exit code for the program is returned as an integer value.
@param argc : The number of parameters being passed in.
@param argv : The alphanumeric parameters from the command line.
@return int : The exit code of the system*/
int main(int argc, char * argv[])
{
	//Load the images from disk
	vector<Mat> templates = loadTemplateImages(vector<Mat>());
	vector<Mat> masks = loadTemplateMasks(vector<Mat>());
	vector<Mat> inputImages = loadInputImages(vector<Mat>());

	//Keeps track of which templates are in which images
	vector< vector<int> > templatesInImage;
	vector<PointVal> resultPointVals;

	//For the number of images in the database
	for (int i = 0; i < inputImages.size(); i++)
	{
		templatesInImage.push_back(vector<int>());

		//For the number of templates
		for (int j = 0; j < templates.size(); j++)
		{
			//initialize an empty PointVal object
			PointVal bestMatch(Point(0, 0), DBL_MIN, Size(0, 0));

			//Number of templates in the image
			string name = "output" + to_string(i) + ".jpg";
			Mat templateParam = templates[j].clone();
			Mat maskParam = masks[j].clone();
			PointVal pv(Point(0,0), 0.0, Size(0, 0));
			templatesInImage = templateInImage(inputImages[i], templateParam, maskParam, bestMatch, templatesInImage, j, name, pv);
			resultPointVals.push_back(pv);
		}
	}

	//Print the results of the template search to disk in a text file
	printResultsToFile(templatesInImage, resultPointVals);

	cin.ignore();

	//Finished Execution
	return 0;
}
//EOF