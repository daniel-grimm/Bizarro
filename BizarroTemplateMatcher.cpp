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
#include <filesystem>
#include <windows.h>
#include <stdio.h>
#include "PointVal.cpp"

//Namespaces
using namespace cv;
using namespace std;

bool debug;

/*Preconditions: An empty vector<Mat>& is passed into the method
for the purpose of containing all of the template images.
Postconditions: The template images are returned as a vector<Mat>&.
@param vector : An empty vector<Mat>& which will contain an array of templates.
@return vector<Mat>& : The array of templates from the database.*/
vector<Mat>& loadTemplateImages(vector<Mat>& vector, int numberOfTemplates)
{
	//For the number of templates in the database
	char character = 'a';

	//Read Image into a Mat object
	for (int i = 0; i < numberOfTemplates; i++)
	{
		//Create the name of the file
		string templateName = ".//templates//";
		templateName.push_back(character);
		templateName += ".jpg";

		//Read in the image
		Mat templateImage = imread(templateName, CV_LOAD_IMAGE_COLOR);

		//add the image to the vector
		vector.push_back(templateImage);

		//increment the name of the template
		character = character++;
	}

	//return the vector of template images
	return vector;
}

/*Preconditions: An empty vector<Mat>& is passed into the method
for the purpose of containing all of the template images.
Postconditions: The template images are returned as a vector<Mat>&.
@param vector : An empty vector<Mat>& which will contain an array of templates.
@return vector<Mat>& : The array of templates from the database.*/
vector<Mat>& loadTemplateMasks(vector<Mat>& vector, int numberOfTemplates)
{
	//For the number of templates (masks) in the database
	char character = 'a';

	//Read Image into a Mat object
	for (int i = 0; i < numberOfTemplates; i++)
	{
		//Create the name of the file
		string templateName = ".//templates//";
		templateName.push_back(character);
		templateName += "_mask.jpg";

		//Read in the image
		Mat templateMaskImage = imread(templateName, IMREAD_GRAYSCALE);

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
vector<Mat>& loadInputImages(vector<Mat>& vector, int numberOfComics)
{
	//For the number of images in the database
	for (int i = 0; i < numberOfComics; i++)
	{
		//Read Image into a Mat object
		string imageName = ".//comics//comic" + to_string(i) + ".jpg";
		Mat image = imread(imageName, CV_LOAD_IMAGE_COLOR);

		//add the image to the vector
		vector.push_back(image);
	}
	return vector;
}

///*Preconditions: A non-null image is passed into the method for scaling up.
//Postconditions: The same image that was passed in is returned but has been scaled up.
//@param templateImage : The image that is going to be scaled.
//@return Mat& : The original image but larger.*/
//Mat& changeTemplateScale(Mat& templateImage, double scalingFactor)
//{
//	//resize the image
//	resize(templateImage, templateImage, Size(), scalingFactor, scalingFactor, INTER_LINEAR);
//
//	//return the scaled image
//	return templateImage;
//}

/*Preconditions: A non-null image is passed into the method for scaling up.
Postconditions: The same image that was passed in is returned but has been scaled up.
@param templateImage : The image that is going to be scaled.
@return Mat& : The original image but larger.*/
Mat changeTemplateScale(const Mat& templateImage, double scalingFactor)
{
	Mat output = templateImage.clone();

	//resize the image
	resize(output, output, Size(), scalingFactor, scalingFactor, INTER_LINEAR);

	//return the scaled image
	return output;
}

/*
Code from:
https://stackoverflow.com/questions/22041699/rotate-an-image-without-cropping-in-opencv-in-c

Preconditions: A non-null image is passed into the method for rotation.
Postconditions: The same image is returned by has been rotated slightly.
@param templateImage : The image being rotated.
@return Mat& : The rotated image.*/
Mat changeTemplateRotation(const Mat templateImage, double degreesOfRotation)
{
	Mat rotated = templateImage.clone();

	//Rotate the template
	Point2f pictureCenter(rotated.rows / 2.0, rotated.cols / 2.0);
	Mat rotationMatrix = getRotationMatrix2D(pictureCenter, degreesOfRotation, 1.0);
	Rect boundingBox = RotatedRect(pictureCenter, rotated.size(), degreesOfRotation).boundingRect();


	// adjust transformation matrix
	rotationMatrix.at<double>(0, 2) += boundingBox.width / 2.0 - pictureCenter.x;
	rotationMatrix.at<double>(1, 2) += boundingBox.height / 2.0 - pictureCenter.y;

	warpAffine(templateImage, rotated, rotationMatrix, boundingBox.size());

	//return the rotated image
	return rotated;
}

void drawOutlinedText(Mat& inputImage, string text, Point point, Scalar back, Scalar fore) {

	cv::putText(inputImage, text, Point(point.x - 1, point.y - 1), FONT_HERSHEY_COMPLEX_SMALL, 1.0, back, 1, CV_AA);
	cv::putText(inputImage, text, Point(point.x - 1, point.y + 1), FONT_HERSHEY_COMPLEX_SMALL, 1.0, back, 1, CV_AA);
	cv::putText(inputImage, text, Point(point.x + 1, point.y - 1), FONT_HERSHEY_COMPLEX_SMALL, 1.0, back, 1, CV_AA);
	cv::putText(inputImage, text, Point(point.x + 1, point.y + 1), FONT_HERSHEY_COMPLEX_SMALL, 1.0, back, 1, CV_AA);
	cv::putText(inputImage, text, point, FONT_HERSHEY_COMPLEX_SMALL, 1.0, fore, 1, CV_AA);

}

/*Preconditions:
Postconditions:
*/
void drawGreenBox(Mat& inputImage, const string label, const PointVal& pv, const string name, const double scalingFactor, const int scalingIterations,
	const double degreesOfRotation, const double thresholdValue, const double maxThreshold, const bool templatePenalty, const double templatePenaltyModifier)
{
	Point point = pv.point;
	Size size = pv.size;

	rectangle(inputImage, point, Point(point.x + size.width, point.y + size.height), Scalar(0, 255, 127), 2);
	drawOutlinedText(inputImage, label, point, Scalar(0, 0, 0), Scalar(0, 255, 127));

	stringstream ss;
	ss << ".\\outputs\\";

	CreateDirectory(ss.str().c_str(), NULL);

	ss << "scalingFactor-" << scalingFactor << "_scalingIterations-" << scalingIterations << "_degreesOfRotation-" << degreesOfRotation << 
		"_thresholdValue-" << thresholdValue << "_maxThreshold-" << maxThreshold << "_templatePenalty-" << templatePenalty <<
		"_templatePenaltyModifier-" << templatePenaltyModifier << "\\";

	CreateDirectory(ss.str().c_str(), NULL);

	ss << name + ".jpg";

	imwrite(ss.str(), inputImage);
}

/*Preconditions:
Postconditions:
*/
void drawGreenBox(const Mat& inputImage, const PointVal& pv, string name, char c)
{
	Point point = pv.point;
	Size size = pv.size;

	Mat img = inputImage.clone();
	cvtColor(img, img, COLOR_GRAY2BGR);

	rectangle(img, point, Point(point.x + size.width, point.y + size.height), Scalar(0, 255, 127), 2);

	stringstream ss;
	ss << ".\\outputs\\";

	CreateDirectory(ss.str().c_str(), NULL);

	ss << "template_" << c << "\\";

	CreateDirectory(ss.str().c_str(), NULL);
		
	ss << name + "_template_" + ".jpg";
	imwrite(ss.str(), img);

	/*stringstream ss;
	ss << "template_" << c << "_" + name + "_" << point.x << "_" << point.y << ".jpg";
	imwrite(ss.str(), img);*/
}

void overlayTemplateImage(Mat src, Mat overlay, const PointVal& pv, string name, char c) {
	
	Size size = pv.size;

	Mat flatGrey(src.rows, src.cols, CV_8UC1, Scalar(128, 128, 128));
	overlay.copyTo(flatGrey(Rect(pv.point.x, pv.point.y, size.width, size.height)));
	
	double alpha = 0.75;
	double beta = 1.0;
	Mat output;
	addWeighted(flatGrey, alpha, src, beta, 0.0, output);

	stringstream ss;
	ss << ".\\outputs\\";

	CreateDirectory(ss.str().c_str(), NULL);

	ss << "template_" << c << "\\";

	CreateDirectory(ss.str().c_str(), NULL);

	ss << name + "_template_" + ".jpg";
	imwrite(ss.str(), output);

}

/*Preconditions: Non-null values for an input image and template image are
passed into the method. The dimensions of the input image are assumed to be
of equal or larger size than the template.
Postconditions: The number of times the template appears in the image is returned.
@param inputImage : The image being searched over.
@param templateImage : The template that is attempting to be found in the image.
@param bestMatch : The current best match for the template in the image.
@return PointVal&: The point and value of the best match of the template in the image.*/
void slideTemplateOverImage(const Mat& inputImage, const Mat& templateImage, const Mat& maskImage, PointVal& bestMatch, const double templateRatio, 
	bool& changed, bool& smallStep, const bool templatePenalty, const double templatePenaltyModifier)
{

	//Create the result matrix
	int resultRows = inputImage.rows - templateImage.rows + 1;
	int resultCols = inputImage.cols - templateImage.cols + 1;
	Mat result(resultRows, resultCols, CV_8UC1);

	//Match the template to the image
	matchTemplate(inputImage, templateImage, result, CV_TM_CCORR_NORMED, maskImage);

	//Get the minimumLocation of the result matrix
	double minimum, maximum;
	Point minLocation, maxLocation;
	minMaxLoc(result, &minimum, &maximum, &minLocation, &maxLocation, Mat());

	//Retrieve the Point and Value of the minimum location
	if (templatePenalty) maximum = maximum * templateRatio;
	PointVal testPoint(maxLocation, maximum, Size(templateImage.cols, templateImage.rows));

	//If a better match is found, update the new point
	if (testPoint > bestMatch && abs(maxLocation.x - bestMatch.point.x) > 4 && abs(maxLocation.y - bestMatch.point.y) > 4)
	{
		if (debug) {
			cout << "updating best match from " << bestMatch.doubleVal << " to " << testPoint.doubleVal << "	-	a " << (testPoint.doubleVal - bestMatch.doubleVal) << " increase" << endl;
			cout << "location from (" << bestMatch.point.x << ", " << bestMatch.point.y << ") to (" << maxLocation.x << ", " << maxLocation.y << ")" << endl;
			cout << "template size: " << templateImage.size << endl;
		}

		if ((testPoint.doubleVal - bestMatch.doubleVal) < 0.01) {
			smallStep = true;
			if (debug) cout << "small step: " << testPoint.doubleVal - bestMatch.doubleVal << endl;
		}
		bestMatch = testPoint;
		changed = true;
	}
}

Mat findEdges(const Mat& image) {

	Mat img = image.clone();

	//create grey image
	cvtColor(img, img, COLOR_BGR2GRAY);

	//blur it up
	GaussianBlur(img, img, Size(5, 5), 2.0, 2.0);

	//Sobel it up
	Sobel(img, img, CV_8U, 1, 0);

	//return the edge image
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
PointVal templateInImage(const Mat& inputImage, Mat& overlayImage, const Mat& templateImage, const Mat& maskImage, vector<string> labels, PointVal& bestMatch,
	vector< vector<int> >& templatesInImage, vector< vector<PointVal> >& pointValuesInImage, const int templateNumber, const string name, const double scalingFactor,
	const int scalingIterations, const double& degreesOfRotation, const double thresholdValue, const double maxThreshold, const bool templatePenalty, 
	const double templatePenaltyModifier)
{
	//create sobel-filtered edge images of the template image
	//and the search image
	Mat image = findEdges(inputImage);

	//initialize an empty pointval
	PointVal test(Point(0, 0), DBL_MIN, Size(0, 0));

	//define constants
	double templateRatio = 1;
	bool changed;
	bool threshold = false;
	bool smallStep = false;
	bool match = false;

	//make an array to hold all the possible rotations of the template image
	vector<Mat> templates;
	vector<Mat> masks;

	vector<Mat> flippedTemplates;
	vector<Mat> flippedMasks;

	//create a mat to hold the scaled version of the template image, so
	//we don't rescale an already scaled version and lose resolution quality
	Mat scaled = templateImage;
	Mat scaledMask = maskImage;

	//for all scales of the template image
	for (int i = 1; i <= scalingIterations; i++) {

		if (threshold) break;

		scaled = findEdges(scaled);

	/*	stringstream ss;
		ss << "template_" << char(templateNumber + 97) << "_" << to_string(i) << ".jpg";
		imwrite(ss.str(), scaled);*/

		//create all the rotations
		for (int j = 0; j < 360 / degreesOfRotation; j++)
		{
			//Rotate the image 360 / degreesOfRotation times and store each in
			//the templates array
			templates.push_back(changeTemplateRotation(scaled, degreesOfRotation * j));
			masks.push_back(changeTemplateRotation(scaledMask, degreesOfRotation * j));

			Mat flipped = templates[j];
			Mat flippedMask = masks[j];
			flip(templates[j], flipped, 1);
			flip(masks[j], flippedMask, 1);

			flippedTemplates.push_back(flipped);
			flippedMasks.push_back(flippedMask);
		}

		//for all rotations of the template image
		for (int j = 0; j < 360 / degreesOfRotation; j++)
		{
			changed = false;
			smallStep = false;

			//Find the best match of the template in the image
			slideTemplateOverImage(image, templates[j], masks[j], bestMatch, templateRatio, changed, smallStep, templatePenalty, templatePenaltyModifier);
			
			if (changed) {

				if (debug) {
					//produce an interim "green box" image to check on the progress of the template matching
					overlayTemplateImage(image, templates[j], bestMatch, name + "_" + to_string(i) + "rotated_by_" + to_string(j * degreesOfRotation) + "_flipped_template", char(templateNumber + 97));
					//drawGreenBox(image, bestMatch, name + "_" + to_string(i) + "rotated_by_" + to_string(j * degreesOfRotation), char(templateNumber + 97));
				}

				if (bestMatch.doubleVal > thresholdValue && smallStep) {
					threshold = true;
					match = true;
					break;
				}
				else if (bestMatch.doubleVal > thresholdValue + maxThreshold) {
					threshold = true;
					break;
				}
			}

			changed = false;
			smallStep = false;

			slideTemplateOverImage(image, flippedTemplates[j], flippedMasks[j], bestMatch, templateRatio, changed, smallStep, templatePenalty, templatePenaltyModifier);

			if (changed) {

				if (debug) {
					//produce an interim "green box" image to check on the progress of the template matching
					overlayTemplateImage(image, flippedTemplates[j], bestMatch, name + "_" + to_string(i) + "rotated_by_" + to_string(j * degreesOfRotation) + "_flipped_template", char(templateNumber + 97));
					//drawGreenBox(image, bestMatch, name + "_" + to_string(i) + "rotated_by_" + to_string(j * degreesOfRotation) + "_flipped_template", char(templateNumber + 97));
				}

				if (bestMatch.doubleVal > thresholdValue && smallStep) {
					threshold = true;
					match = true;
					break;
				}
				else if (bestMatch.doubleVal > thresholdValue + maxThreshold) {
					threshold = true;
					break;
				}
			}
		}

		//rescale the template
		double newScale = pow(scalingFactor, i);
		scaled = changeTemplateScale(templateImage, newScale);
		scaledMask = changeTemplateScale(maskImage, newScale);
		templateRatio = newScale + (i * templatePenaltyModifier);

		templates.clear();
		masks.clear();
		flippedTemplates.clear();
		flippedMasks.clear();

	}

	PointVal retVal(Point(0,0),0,Size());

	//draw a green box over the best match of the template
	if (match)
	{

		drawGreenBox(overlayImage, labels[templateNumber], bestMatch, name, scalingFactor, scalingIterations, degreesOfRotation,
			thresholdValue, maxThreshold, templatePenalty, templatePenaltyModifier);

		//Add the best match to the vector
		templatesInImage[templatesInImage.size() - 1].push_back(templateNumber);
		pointValuesInImage[templatesInImage.size() - 1].push_back(bestMatch);

		if (debug) {
			//print the template's best match's location and value
			cout << "MATCH: Point: " << bestMatch.point << " Value: " << bestMatch.doubleVal << endl << endl;
		}
	}
	else {

		if (debug) {
			cout << "NO MATCH" << endl;
		}

	}
	
	

	return retVal;
}

/*Preconditions: An empty vector<string> is passed into the method which
will contain the name of all the templates
Postconditions: The names of all the template images are loaded into the vector
at each index.
@param templateNames : An empty vector of strings.
@return vector<string>& : The list of template names.*/
vector<string>& loadStrings(vector<string>& templateNames, string filename)
{
	//open a file to read in the template names
	string title;
	ifstream file(filename);

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
void printResultsToFile(vector< vector<int> >& vectorOfMatchedTemplates, vector< vector<PointVal> > pointValuesInImage)
{
	vector<string> templateNames;
	templateNames = loadStrings(templateNames, "templateNames.txt");

	ofstream out("results.txt");

	for (int i = 0; i < vectorOfMatchedTemplates.size(); i++)
	{
		out << "Image comic" << to_string(i) << ".jpg:" << endl;

		if (vectorOfMatchedTemplates[i].empty()) {
			out << "No templates found." << endl;
		}

		for (int j = 0; j < vectorOfMatchedTemplates[i].size(); j++)
		{
			out << templateNames[vectorOfMatchedTemplates[i][j]] << " - matched at: " << pointValuesInImage[i][j].point << " with certainty: " << pointValuesInImage[i][j].doubleVal << endl;
		}
		out << endl;
	}

	out.close();
}

/*Preconditions: This program must be run from the command line.
Postconditions: The exit code for the program is returned as an integer value.
@param argc : The number of parameters being passed in.
@param argv : The alphanumeric parameters from the command line.
@return int : The exit code of the system*/
int main(int argc, char * argv[])
{

	const int numberOfTemplates = atoi(argv[1]);
	const int numberOfComics = atoi(argv[2]);
	const double scalingFactor = atof(argv[3]);
	const int scalingIterations = atoi(argv[4]);
	const double degreesOfRotation = atof(argv[5]);
	const double thresholdValue = atof(argv[6]);
	const double maxThreshold = atof(argv[7]);
	string templatePenaltyString = argv[8];
	const bool templatePenalty = (templatePenaltyString == "true");
	const double templatePenaltyModifier = atof(argv[9]);
	string debugString = argv[10];
	debug = (debugString == "true");

	//Load the images from disk
	vector<Mat> templates = loadTemplateImages(vector<Mat>(), numberOfTemplates);
	vector<Mat> masks = loadTemplateMasks(vector<Mat>(), numberOfTemplates);
	vector<Mat> inputImages = loadInputImages(vector<Mat>(), numberOfComics);

	//Keeps track of which templates are in which images
	vector< vector<int> > templatesInImage;
	vector< vector<PointVal> > pointValuesInImage;

	//load in the labels
	vector<string> labels;
	labels = loadStrings(labels, "shortTemplateNames.txt");

	//For the number of images in the database
	for (int i = 0; i < inputImages.size(); i++)
	{
		templatesInImage.push_back(vector<int>());
		pointValuesInImage.push_back(vector<PointVal>());
		Mat output = inputImages[i].clone();

		//For the number of templates
		for (int j = 0; j < templates.size(); j++)
		{
			//initialize an empty PointVal object
			PointVal bestMatch(Point(0, 0), DBL_MIN, Size(0, 0));

			//Number of templates in the image
			string name = "output" + to_string(i);

			Mat templateParam = templates[j].clone();
			Mat maskParam = masks[j].clone();
			

			PointVal match = templateInImage(inputImages[i], output, templateParam, maskParam, labels, bestMatch, templatesInImage, pointValuesInImage, j, name,
				scalingFactor, scalingIterations, degreesOfRotation, thresholdValue, maxThreshold, templatePenalty, templatePenaltyModifier);
		}
	}

	//Print the results of the template search to disk in a text file
	printResultsToFile(templatesInImage, pointValuesInImage);

	cout << endl << "Press Enter to Exit the Program" << endl;
	cin.ignore();

	//Finished Execution
	return 0;
}
//EOF