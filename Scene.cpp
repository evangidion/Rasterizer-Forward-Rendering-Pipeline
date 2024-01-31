#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

void Scene::assignColorToPixel(int i, int j, Color c)
{
	this->image[i][j].r = c.r;
	this->image[i][j].g = c.g;
	this->image[i][j].b = c.b;
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			this->depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				assignColorToPixel(i, j, this->backgroundColor);
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
				this->depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
*/
void Scene::convertPPMToPNG(string ppmFileName)
{
	string command;

	// TODO: Change implementation if necessary.
	// command = "./magick convert " + ppmFileName + " " + ppmFileName + ".png";
	command = "convert " + ppmFileName + " " + ppmFileName + ".png";
	system(command.c_str());
}

/*
	Transformations, clipping, culling, rasterization are done here.
*/

Matrix4 Scene::modelingTransformationMatrix(Mesh* mesh) {
	Matrix4 modelingMatrix = getIdentityMatrix();
	for(int i = 0; i < mesh->numberOfTransformations; i++) {
		if(mesh->transformationTypes[i] == 't') {
			Matrix4 translationMatrix;
			translationMatrix.values[0][0] = 1.0;
			translationMatrix.values[0][1] = 0.0;
			translationMatrix.values[0][2] = 0.0;
			translationMatrix.values[0][3] = translations[mesh->transformationIds[i] - 1]->tx;
			translationMatrix.values[1][0] = 0.0;
			translationMatrix.values[1][1] = 1.0;
			translationMatrix.values[1][2] = 0.0;
			translationMatrix.values[1][3] = translations[mesh->transformationIds[i] - 1]->ty;
			translationMatrix.values[2][0] = 0.0;
			translationMatrix.values[2][1] = 0.0;
			translationMatrix.values[2][2] = 1.0;
			translationMatrix.values[2][3] = translations[mesh->transformationIds[i] - 1]->tz;
			translationMatrix.values[3][0] = 0.0;
			translationMatrix.values[3][1] = 0.0;
			translationMatrix.values[3][2] = 0.0;
			translationMatrix.values[3][3] = 1.0;
			modelingMatrix = multiplyMatrixWithMatrix(translationMatrix, modelingMatrix);
		}
		else if(mesh->transformationTypes[i] == 's') {
			Matrix4 scalingMatrix;
			scalingMatrix.values[0][0] = scalings[mesh->transformationIds[i] - 1]->sx;
			scalingMatrix.values[0][1] = 0.0;
			scalingMatrix.values[0][2] = 0.0;
			scalingMatrix.values[0][3] = 0.0;
			scalingMatrix.values[1][0] = 0.0;
			scalingMatrix.values[1][1] = scalings[mesh->transformationIds[i] - 1]->sy;
			scalingMatrix.values[1][2] = 0.0;
			scalingMatrix.values[1][3] = 0.0;
			scalingMatrix.values[2][0] = 0.0;
			scalingMatrix.values[2][1] = 0.0;
			scalingMatrix.values[2][2] = scalings[mesh->transformationIds[i] - 1]->sz;
			scalingMatrix.values[2][3] = 0.0;
			scalingMatrix.values[3][0] = 0.0;
			scalingMatrix.values[3][1] = 0.0;
			scalingMatrix.values[3][2] = 0.0;
			scalingMatrix.values[3][3] = 1.0;
			modelingMatrix = multiplyMatrixWithMatrix(scalingMatrix, modelingMatrix);
		}
		else if(mesh->transformationTypes[i] == 'r') {
			Vec3 u, v, w;
			u.x = rotations[mesh->transformationIds[i] - 1]->ux;
			u.y = rotations[mesh->transformationIds[i] - 1]->uy;
			u.z = rotations[mesh->transformationIds[i] - 1]->uz;
			u = normalizeVec3(u);
			double minComponent = min(min(ABS(u.x), ABS(u.y)), ABS(u.z));
			if(minComponent == ABS(u.x)) {
				v.x = 0.0;
				v.y = -1.0 * u.z;
				v.z = u.y;
			}
			else if(minComponent == ABS(u.y)) {
				v.x = -1.0 * u.z;
				v.y = 0.0;
				v.z = u.x;
			}
			else {
				v.x = -1.0 * u.y;
				v.y = u.x;
				v.z = 0.0;
			}
			v = normalizeVec3(v);
			w = crossProductVec3(u, v);
			w = normalizeVec3(w);
			Matrix4 M, MT;
			M.values[0][0] = u.x;
			M.values[0][1] = u.y;
			M.values[0][2] = u.z;
			M.values[0][3] = 0.0;
			M.values[1][0] = v.x;
			M.values[1][1] = v.y;
			M.values[1][2] = v.z;
			M.values[1][3] = 0.0;
			M.values[2][0] = w.x;
			M.values[2][1] = w.y;
			M.values[2][2] = w.z;
			M.values[2][3] = 0.0;
			M.values[3][0] = 0.0;
			M.values[3][1] = 0.0;
			M.values[3][2] = 0.0;
			M.values[3][3] = 1.0;
			MT.values[0][0] = u.x;
			MT.values[0][1] = v.x;
			MT.values[0][2] = w.x;
			MT.values[0][3] = 0.0;
			MT.values[1][0] = u.y;
			MT.values[1][1] = v.y;
			MT.values[1][2] = w.y;
			MT.values[1][3] = 0.0;
			MT.values[2][0] = u.z;
			MT.values[2][1] = v.z;
			MT.values[2][2] = w.z;
			MT.values[2][3] = 0.0;
			MT.values[3][0] = 0.0;
			MT.values[3][1] = 0.0;
			MT.values[3][2] = 0.0;
			MT.values[3][3] = 1.0;
			double radian = (rotations[mesh->transformationIds[i] - 1]->angle * M_PI) / 180.0;
			Matrix4 rotationMatrix;
			rotationMatrix.values[0][0] = 1.0;
			rotationMatrix.values[0][1] = 0.0; 
			rotationMatrix.values[0][2] = 0.0;
			rotationMatrix.values[0][3] = 0.0;
			rotationMatrix.values[1][0] = 0.0;
			rotationMatrix.values[1][1] = cos(radian);
			rotationMatrix.values[1][2] = -1.0 * sin(radian);
			rotationMatrix.values[1][3] = 0.0;
			rotationMatrix.values[2][0] = 0.0;
			rotationMatrix.values[2][1] = sin(radian);
			rotationMatrix.values[2][2] = cos(radian);
			rotationMatrix.values[2][3] = 0.0;
			rotationMatrix.values[3][0] = 0.0;
			rotationMatrix.values[3][1] = 0.0;
			rotationMatrix.values[3][2] = 0.0;
			rotationMatrix.values[3][3] = 1.0;
			rotationMatrix = multiplyMatrixWithMatrix(rotationMatrix, M);
			rotationMatrix = multiplyMatrixWithMatrix(MT, rotationMatrix);
			modelingMatrix = multiplyMatrixWithMatrix(rotationMatrix, modelingMatrix);
		}
	}
	return modelingMatrix;
}

Matrix4 cameraTransformationMatrix(Camera *camera) {
	Matrix4 cameraMatrix;
	cameraMatrix.values[0][0] = camera->u.x;
	cameraMatrix.values[0][1] = camera->u.y;
	cameraMatrix.values[0][2] = camera->u.z;
	cameraMatrix.values[0][3] = -1.0 * (camera->u.x * camera->position.x + camera->u.y * camera->position.y + camera->u.z * camera->position.z); 
	cameraMatrix.values[1][0] = camera->v.x;
	cameraMatrix.values[1][1] = camera->v.y;
	cameraMatrix.values[1][2] = camera->v.z;
	cameraMatrix.values[1][3] = -1.0 * (camera->v.x * camera->position.x + camera->v.y * camera->position.y + camera->v.z * camera->position.z);
	cameraMatrix.values[2][0] = camera->w.x;
	cameraMatrix.values[2][1] = camera->w.y;
	cameraMatrix.values[2][2] = camera->w.z;
	cameraMatrix.values[2][3] = -1.0 * (camera->w.x * camera->position.x + camera->w.y * camera->position.y + camera->w.z * camera->position.z);
	cameraMatrix.values[3][0] = 0.0;
	cameraMatrix.values[3][1] = 0.0;
	cameraMatrix.values[3][2] = 0.0;
	cameraMatrix.values[3][3] = 1.0;
	return cameraMatrix;
}

Matrix4 projectionTransformationMatrix(Camera *camera) {
	if(camera->projectionType) {
		Matrix4 perspectiveMatrix;
		perspectiveMatrix.values[0][0] = (2.0 * camera->near) / (camera->right - camera->left);
		perspectiveMatrix.values[0][1] = 0.0;
		perspectiveMatrix.values[0][2] = (camera->right + camera->left) / (camera->right - camera->left);
		perspectiveMatrix.values[0][3] = 0.0;
		perspectiveMatrix.values[1][0] = 0.0;
		perspectiveMatrix.values[1][1] = (2.0 * camera->near) / (camera->top - camera->bottom);
		perspectiveMatrix.values[1][2] = (camera->top + camera->bottom) / (camera->top - camera->bottom);
		perspectiveMatrix.values[1][3] = 0.0;
		perspectiveMatrix.values[2][0] = 0.0;
		perspectiveMatrix.values[2][1] = 0.0;
		perspectiveMatrix.values[2][2] = (camera->far + camera->near) / (camera->near - camera->far);
		perspectiveMatrix.values[2][3] = (2.0 * camera->far * camera->near) / (camera->near - camera->far);
		perspectiveMatrix.values[3][0] = 0.0;
		perspectiveMatrix.values[3][1] = 0.0;
		perspectiveMatrix.values[3][2] = -1.0;
		perspectiveMatrix.values[3][3] = 0.0;
		return perspectiveMatrix;
	}
	else {
		Matrix4 orthographicMatrix;
		orthographicMatrix.values[0][0] = 2.0 / (camera->right - camera->left);
		orthographicMatrix.values[0][1] = 0.0;
		orthographicMatrix.values[0][2] = 0.0;
		orthographicMatrix.values[0][3] = (camera->right + camera->left) / (camera->left - camera->right);
		orthographicMatrix.values[1][0] = 0.0;
		orthographicMatrix.values[1][1] = 2.0 / (camera->top - camera->bottom);
		orthographicMatrix.values[1][2] = 0.0;
		orthographicMatrix.values[1][3] = (camera->top + camera->bottom) / (camera->bottom - camera->top);
		orthographicMatrix.values[2][0] = 0.0;
		orthographicMatrix.values[2][1] = 0.0;
		orthographicMatrix.values[2][2] = 2.0 / (camera->near - camera->far);
		orthographicMatrix.values[2][3] = (camera->far + camera->near) / (camera->near - camera->far);
		orthographicMatrix.values[3][0] = 0.0;
		orthographicMatrix.values[3][1] = 0.0;
		orthographicMatrix.values[3][2] = 0.0;
		orthographicMatrix.values[3][3] = 1.0;
		return orthographicMatrix;
	}
}

void perspectiveDivide(Vec4 &p) {
	p.x /= p.t;
	p.y /= p.t;
	p.z /= p.t;
	p.t = 1.0;
}

Matrix4 viewportTransformationMatrix(Camera *camera) {
	Matrix4 viewportMatrix;
	viewportMatrix.values[0][0] = camera->horRes / 2.0;
	viewportMatrix.values[0][1] = 0.0;
	viewportMatrix.values[0][2] = 0.0;
	viewportMatrix.values[0][3] = (camera->horRes - 1.0) / 2.0;
	viewportMatrix.values[1][0] = 0.0;
	viewportMatrix.values[1][1] = camera->verRes / 2.0;
	viewportMatrix.values[1][2] = 0.0;
	viewportMatrix.values[1][3] = (camera->verRes - 1.0) / 2.0;
	viewportMatrix.values[2][0] = 0.0;
	viewportMatrix.values[2][1] = 0.0;
	viewportMatrix.values[2][2] = 0.5;
	viewportMatrix.values[2][3] = 0.5;
	viewportMatrix.values[3][0] = 0.0;
	viewportMatrix.values[3][1] = 0.0;
	viewportMatrix.values[3][2] = 0.0;
	viewportMatrix.values[3][3] = 0.0;
	return viewportMatrix;
}

bool visible(double den, double num, double &tE, double &tL) {
	if(den > 0) {
		double t = num / den;
		if(t > tL) return false;
		if(t > tE) tE = t;
	}
	else if(den < 0) {
		double t = num / den;
		if(t < tE) return false;
		if(t < tL) tL = t;
	}
	else if(num > 0) return false;
	return true;
}

bool liangBarskyClipping(Vec4 &p1, Color &c1, Vec4 &p2, Color &c2) {
	double tE = 0.0;
	double tL = 1.0;
	bool vision = false;
	double dx = p2.x - p1.x;
	double dy = p2.y - p1.y;
	double dz = p2.z - p1.z;
	Color dc;
	dc.r = c2.r - c1.r;
	dc.g = c2.g - c1.g;
	dc.b = c2.b - c1.b;
	double w = ABS(p1.t) >= ABS(p2.t) ? ABS(p1.t) : ABS(p2.t);
	if(visible(dx, -w - p1.x, tE, tL)) {
		if(visible(-dx, p1.x - w, tE, tL)) {
			if(visible(dy, -w - p1.y, tE, tL)) {
				if(visible(-dy, p1.y - w, tE, tL)) {
					if(visible(dz, -w - p1.z, tE, tL)) {
						if(visible(-dz, p1.z - w, tE, tL)) {
							vision = true;
							if(tL < 1.0) {
								p2.x = p1.x + dx * tL;
								p2.y = p1.y + dy * tL;
								p2.z = p1.z + dz * tL;
								c2.r = c1.r + dc.r * tL;
								c2.g = c1.g + dc.g * tL;
								c2.b = c1.b + dc.b * tL;
							}
							if(tE > 0.0) {
								p1.x += dx * tE;
								p1.y += dy * tE;
								p1.z += dz * tE;
								c1.r += dc.r * tE;
								c1.g += dc.g * tE;
								c1.b += dc.b * tE;
							}
						}
					}
				}
			}
		}
	}
	return vision;
}

void midpointAlgorithm(vector<vector<Color>> &image, vector<vector<double>> &depth, Vec4 p0, Color c0, Vec4 p1, Color c1) {
	if(ABS(p1.y - p0.y) <= ABS(p1.x - p0.x)) {
		if(p1.x - p0.x >= 0) {
			int dx = p1.x - p0.x;
			int dy = p1.y - p0.y;
			if(dy >= 0) {
				int y = p0.y;
				int d = 2 * (-dy) + dx;
				Color c = c0;
				Color dc;
				dc.r = (c1.r - c0.r) / dx;
				dc.g = (c1.g - c0.g) / dx;
				dc.b = (c1.b - c0.b) / dx;
				double z = p0.z;
				double dz = (p1.z - p0.z) / dx;
				for(int x = p0.x; x <= (int) p1.x; x++) {
					if(0 <= x && x < image.size() && 0 <= y && y < image[0].size() && z <= depth[x][y]) {
						image[x][y] = c;
						depth[x][y] = z;
					}
					if(d < 0) {
						y++;
						d += 2 * (-dy + dx);
					}
					else d += 2 * (-dy);
					c.r += dc.r;
					c.g += dc.g;
					c.b += dc.b;	
					z += dz;			
				}				
			}
			else {
				int y = p0.y;
				int d = 2 * (-dy) - dx;
				Color c = c0;
				Color dc;
				dc.r = (c1.r - c0.r) / dx;
				dc.g = (c1.g - c0.g) / dx;
				dc.b = (c1.b - c0.b) / dx;
				double z = p0.z;
				double dz = (p1.z - p0.z) / dx;
				for(int x = p0.x; x <= (int) p1.x; x++) {
					if(0 <= x && x < image.size() && 0 <= y && y < image[0].size() && z <= depth[x][y]) {
						image[x][y] = c;
						depth[x][y] = z;
					}
					if(d > 0) {
						y--;
						d += 2 * (-dy - dx);
					}
					else d += 2 * (-dy);
					c.r += dc.r;
					c.g += dc.g;
					c.b += dc.b;	
					z += dz;			
				}
			}
		}
		else {
			Vec4 p = p0;
			p0 = p1;
			p1 = p;
			Color temp = c0;
			c0 = c1;
			c1 = temp;
			int dx = p1.x - p0.x;
			int dy = p1.y - p0.y;
			if(dy >= 0) {
				int y = p0.y;
				int d = 2 * (-dy) + dx;
				Color c = c0;
				Color dc;
				dc.r = (c1.r - c0.r) / dx;
				dc.g = (c1.g - c0.g) / dx;
				dc.b = (c1.b - c0.b) / dx;
				double z = p0.z;
				double dz = (p1.z - p0.z) / dx;
				for(int x = p0.x; x <= (int) p1.x; x++) {
					if(0 <= x && x < image.size() && 0 <= y && y < image[0].size() && z <= depth[x][y]) {
						image[x][y] = c;
						depth[x][y] = z;
					}
					if(d < 0) {
						y++;
						d += 2 * (-dy + dx);
					}
					else d += 2 * (-dy);
					c.r += dc.r;
					c.g += dc.g;
					c.b += dc.b;	
					z += dz;			
				}				
			}
			else {
				int y = p0.y;
				int d = 2 * (-dy) - dx;
				Color c = c0;
				Color dc;
				dc.r = (c1.r - c0.r) / dx;
				dc.g = (c1.g - c0.g) / dx;
				dc.b = (c1.b - c0.b) / dx;
				double z = p0.z;
				double dz = (p1.z - p0.z) / dx;
				for(int x = p0.x; x <= (int) p1.x; x++) {
					if(0 <= x && x < image.size() && 0 <= y && y < image[0].size() && z <= depth[x][y]) {
						image[x][y] = c;
						depth[x][y] = z;
					}
					if(d > 0) {
						y--;
						d += 2 * (-dy - dx);
					}
					else d += 2 * (-dy);
					c.r += dc.r;
					c.g += dc.g;
					c.b += dc.b;	
					z += dz;			
				}
			}			
		}
	}
	else {
		if(p1.y - p0.y >= 0) {
			int dx = p1.x - p0.x;
			int dy = p1.y - p0.y;
			if(dx >= 0) {
				int x = p0.x;
				int d = 2 * (-dx) + dy;
				Color c = c0;
				Color dc;
				dc.r = (c1.r - c0.r) / dy;
				dc.g = (c1.g - c0.g) / dy;
				dc.b = (c1.b - c0.b) / dy;
				double z = p0.z;
				double dz = (p1.z - p0.z) / dy;
				for(int y = p0.y; y <= (int) p1.y; y++) {
					if(0 <= x && x < image.size() && 0 <= y && y < image[0].size() && z <= depth[x][y]) {
						image[x][y] = c;
						depth[x][y] = z;
					}
					if(d < 0) {
						x++;
						d += 2 * (-dx + dy);
					}
					else d += 2 * (-dx);
					c.r += dc.r;
					c.g += dc.g;
					c.b += dc.b;	
					z += dz;			
				}				
			}
			else {
				int x = p0.x;
				int d = 2 * (-dx) - dy;
				Color c = c0;
				Color dc;
				dc.r = (c1.r - c0.r) / dy;
				dc.g = (c1.g - c0.g) / dy;
				dc.b = (c1.b - c0.b) / dy;
				double z = p0.z;
				double dz = (p1.z - p0.z) / dy;
				for(int y = p0.y; y <= (int) p1.y; y++) {
					if(0 <= x && x < image.size() && 0 <= y && y < image[0].size() && z <= depth[x][y]) {
						image[x][y] = c;
						depth[x][y] = z;
					}
					if(d > 0) {
						x--;
						d += 2 * (-dx - dy);
					}
					else d += 2 * (-dx);
					c.r += dc.r;
					c.g += dc.g;
					c.b += dc.b;
					z += dz;				
				}
			}
		}
		else {
			Vec4 p = p0;
			p0 = p1;
			p1 = p;
			Color temp = c0;
			c0 = c1;
			c1 = temp;
			int dx = p1.x - p0.x;
			int dy = p1.y - p0.y;
			if(dx >= 0) {
				int x = p0.x;
				int d = 2 * (-dx) + dy;
				Color c = c0;
				Color dc;
				dc.r = (c1.r - c0.r) / dy;
				dc.g = (c1.g - c0.g) / dy;
				dc.b = (c1.b - c0.b) / dy;
				double z = p0.z;
				double dz = (p1.z - p0.z) / dy;
				for(int y = p0.y; y <= (int) p1.y; y++) {
					if(0 <= x && x < image.size() && 0 <= y && y < image[0].size() && z <= depth[x][y]) {
						image[x][y] = c;
						depth[x][y] = z;
					}
					if(d < 0) {
						x++;
						d += 2 * (-dx + dy);
					}
					else d += 2 * (-dx);
					c.r += dc.r;
					c.g += dc.g;
					c.b += dc.b;	
					z += dz;			
				}				
			}
			else {
				int x = p0.x;
				int d = 2 * (-dx) - dy;
				Color c = c0;
				Color dc;
				dc.r = (c1.r - c0.r) / dy;
				dc.g = (c1.g - c0.g) / dy;
				dc.b = (c1.b - c0.b) / dy;
				double z = p0.z;
				double dz = (p1.z - p0.z) / dy;
				for(int y = p0.y; y <= (int) p1.y; y++) {
					if(0 <= x && x < image.size() && 0 <= y && y < image[0].size() && z <= depth[x][y]) {
						image[x][y] = c;
						depth[x][y] = z;
					}
					if(d > 0) {
						x--;
						d += 2 * (-dx - dy);
					}
					else d += 2 * (-dx);
					c.r += dc.r;
					c.g += dc.g;
					c.b += dc.b;	
					z += dz;			
				}
			}
		}
	}
}

double implicitLineEquation(int x, int y, int x1, int y1, int x2, int y2) {
	return x * (y1 - y2) + y * (x2 - x1) + x1 * y2 - y1 * x2;
}

void triangleRasterization(vector<vector<Color>> &image, vector<vector<double>> &depth, Vec4 p1, Color c1, Vec4 p2, Color c2, Vec4 p3, Color c3) {
	int xMin = p1.x < p2.x ? (int) p1.x : (int) p2.x;
	xMin = xMin < (int) p3.x ? xMin : (int) p3.x;
	int xMax = p1.x > p2.x ? (int) p1.x : (int) p2.x;
	xMax = xMax > (int) p3.x ? xMax : (int) p3.x;
	int yMin = p1.y < p2.y ? (int) p1.y : (int) p2.y;
	yMin = yMin < (int) p3.y ? yMin : (int) p3.y;
	int yMax = p1.y > p2.y ? (int) p1.y : (int) p2.y;
	yMax = yMax > (int) p3.y ? yMax : (int) p3.y;
	double alpha, beta, gama;
	Color color;
	double implicitLineEquation1 = implicitLineEquation(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
	double implicitLineEquation2 = implicitLineEquation(p2.x, p2.y, p3.x, p3.y, p1.x, p1.y);
	double implicitLineEquation3 = implicitLineEquation(p3.x, p3.y, p1.x, p1.y, p2.x, p2.y);
	for(int y = yMin; y <= yMax; y++) {
		for(int x = xMin; x <= xMax; x++) {
			alpha = implicitLineEquation(x, y, p2.x, p2.y, p3.x, p3.y) / implicitLineEquation1;
			beta = implicitLineEquation(x, y, p3.x, p3.y, p1.x, p1.y) / implicitLineEquation2;
			gama = implicitLineEquation(x, y, p1.x, p1.y, p2.x, p2.y) / implicitLineEquation3;
			if(alpha >= 0.0 && beta >= 0.0 && gama >= 0.0) {
				color.r = c1.r * alpha + c2.r * beta + c3.r * gama;
				color.g = c1.g * alpha + c2.g * beta + c3.g * gama;
				color.b = c1.b * alpha + c2.b * beta + c3.b * gama;
				double z = p1.z * alpha + p2.z * beta + p3.z * gama;
				if (0 <= x && x < image.size() && 0 <= y && y < image[0].size() && z <= depth[x][y]) {
					image[x][y] = color;
					depth[x][y] = z;
				}
			}
		}
	}	
}

bool backfaceCulling(Vec3 v1, Vec3 v2, Vec3 v3) {
	Vec3 edge1 = subtractVec3(v2, v1);
	Vec3 edge2 = subtractVec3(v3, v1);
	return dotProductVec3(crossProductVec3(edge2, edge1), v1) > 0;
}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	Matrix4 cameraMatrix = cameraTransformationMatrix(camera);
	Matrix4 projectionMatrix = projectionTransformationMatrix(camera);
	Matrix4 projectionXcameraMatrix = multiplyMatrixWithMatrix(projectionMatrix, cameraMatrix);
	Matrix4 viewportMatrix = viewportTransformationMatrix(camera);
	for(int i = 0; i < meshes.size(); i++) {
		Matrix4 modelingMatrix = modelingTransformationMatrix(meshes[i]);
		Matrix4 projectionXcameraXmodelingMatrix = multiplyMatrixWithMatrix(projectionXcameraMatrix, modelingMatrix);
		for(int j = 0; j < meshes[i]->numberOfTriangles; j++) {
			Vec3 v1 = *vertices[meshes[i]->triangles[j].vertexIds[0] - 1];
			Vec3 v2 = *vertices[meshes[i]->triangles[j].vertexIds[1] - 1];
			Vec3 v3 = *vertices[meshes[i]->triangles[j].vertexIds[2] - 1];
			Vec4 p1, p2, p3;
			p1.x = v1.x;
			p1.y = v1.y;
			p1.z = v1.z;
			p1.t = 1.0;
			p1.colorId = v1.colorId;
			p2.x = v2.x;
			p2.y = v2.y;
			p2.z = v2.z;
			p2.t = 1.0;
			p2.colorId = v2.colorId;
			p3.x = v3.x;
			p3.y = v3.y;
			p3.z = v3.z;
			p3.t = 1.0;
			p3.colorId = v3.colorId;
			p1 = multiplyMatrixWithVec4(projectionXcameraXmodelingMatrix, p1);
			p2 = multiplyMatrixWithVec4(projectionXcameraXmodelingMatrix, p2);
			p3 = multiplyMatrixWithVec4(projectionXcameraXmodelingMatrix, p3);
			v1.x = p1.x;
			v1.y = p1.y;
			v1.z = p1.z;
			v2.x = p2.x;
			v2.y = p2.y;
			v2.z = p2.z;
			v3.x = p3.x;
			v3.y = p3.y;
			v3.z = p3.z;
			if(cullingEnabled && backfaceCulling(v1, v2, v3)) continue;
			if(meshes[i]->type) {
				perspectiveDivide(p1);
				perspectiveDivide(p2);
				perspectiveDivide(p3);
				p1 = multiplyMatrixWithVec4(viewportMatrix, p1);
				p2 = multiplyMatrixWithVec4(viewportMatrix, p2);
				p3 = multiplyMatrixWithVec4(viewportMatrix, p3);
				Color c1 = *colorsOfVertices[p1.colorId - 1];
				Color c2 = *colorsOfVertices[p2.colorId - 1];
				Color c3 = *colorsOfVertices[p3.colorId - 1];
				triangleRasterization(image, depth, p1, c1, p2, c2, p3, c3);
			}
			else {
				Vec4 p1c = p1;
				Vec4 p2c = p2;
				Vec4 p3c = p3;
				Color c1 = *colorsOfVertices[p1.colorId - 1];
				Color c2 = *colorsOfVertices[p2.colorId - 1];
				Color c3 = *colorsOfVertices[p3.colorId - 1];
				Color c1c = c1;
				Color c2c = c2;
				Color c3c = c3;
				if(liangBarskyClipping(p1, c1, p2, c2)) {
					perspectiveDivide(p1);
					perspectiveDivide(p2);
					p1 = multiplyMatrixWithVec4(viewportMatrix, p1);
					p2 = multiplyMatrixWithVec4(viewportMatrix, p2);
					midpointAlgorithm(image, depth, p1, c1, p2, c2);
				}
				if(liangBarskyClipping(p2c, c2c, p3, c3)) {
					perspectiveDivide(p2c);
					perspectiveDivide(p3);
					p2c = multiplyMatrixWithVec4(viewportMatrix, p2c);
					p3 = multiplyMatrixWithVec4(viewportMatrix, p3);
					midpointAlgorithm(image, depth, p2c, c2c, p3, c3);
				}
				if(liangBarskyClipping(p3c, c3c, p1c, c1c)) {
					perspectiveDivide(p3c);
					perspectiveDivide(p1c);
					p3c = multiplyMatrixWithVec4(viewportMatrix, p3c);
					p1c = multiplyMatrixWithVec4(viewportMatrix, p1c);
					midpointAlgorithm(image, depth, p3c, c3c, p1c, c1c);
				}
			}
		}
	}
}
