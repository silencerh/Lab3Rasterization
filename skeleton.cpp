#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <math.h>

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;
using glm::vec2;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
vec3 cameraPos(0, 0, -3.001);
vec3 currentColor;
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
struct Pixel
{
	int x;
	int y;
	float zinv;
	vec3 illumination;
};
struct Vertex
{
	vec3 position;
	vec3 normal;
	vec3 reflectance;
};
vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 1.1f*vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);
//vec3 currentNormal;
//vec3 currentReflectance;
// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();
void VertexShader(const vec3& v, ivec2& p);
void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color);
void DrawPolygonEdges(const vector<vec3>& vertices);
void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels,
	vector<ivec2>& rightPixels);
void DrawRows(const vector<ivec2>& leftPixels, const vector<ivec2>& rightPixels);
void DrawPolygon(const vector<Vertex>& vertices);

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result);
void ComputePolygonRows(
	const vector<Pixel>& vertexPixels,
	vector<Pixel>& leftPixels,
	vector<Pixel>& rightPixels
	);

void DrawPolygonRows(
	const vector<Pixel>& leftPixels,
	const vector<Pixel>& rightPixels
	);
void VertexShader(const Vertex& v, Pixel& p);
void PixelShader(const Pixel& p);

int main(int argc, char* argv[])
{
	LoadTestModel(triangles);
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.

	while (NoQuitMessageSDL())
	{
		Update();
		Draw();
	}

	SDL_SaveBMP(screen, "screenshot.bmp");
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	/*Uint8* keystate = SDL_GetKeyState(0);

	if (keystate[SDLK_UP])
	;

	if( keystate[SDLK_DOWN] )
	;

	if( keystate[SDLK_RIGHT] )
	;

	if( keystate[SDLK_LEFT] )
	;

	if( keystate[SDLK_RSHIFT] )
	;

	if( keystate[SDLK_RCTRL] )
	;

	if( keystate[SDLK_w] )
	;

	if( keystate[SDLK_s] )
	;

	if( keystate[SDLK_d] )
	;

	if( keystate[SDLK_a] )
	;

	if( keystate[SDLK_e] )
	;

	if( keystate[SDLK_q] )
	;*/
}

void Draw()
{
	SDL_FillRect(screen, 0, 0);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for (int y = 0; y < SCREEN_HEIGHT; ++y)
	{
		for (int x = 0; x < SCREEN_WIDTH; ++x)
		{
			depthBuffer[y][x] = 0;
		}
	}
	for (int i = 0; i<triangles.size(); ++i)
	{
		vector<Vertex> vertices(3);

		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;
		//currentNormal = triangles[i].normal;
		//currentReflectance = triangles[i].color;
		for (int j = 0; j < 3; j++)
		{
			vertices[j].normal = triangles[i].normal;
			vertices[j].reflectance = triangles[i].color;
		}
		//part 4
		
		//part 3
		currentColor = triangles[i].color;
		DrawPolygon(vertices);

		//part 2
		//vector<ivec2> projPos(3);
		//vec3 color(1, 1, 1);
		//DrawPolygonEdges(vertices);

		//part 1
		/*for (int v = 0; v < 3; ++v)
		{
			VertexShader(vertices[v], projPos[v]);
			PutPixelSDL(screen, projPos[v].x, projPos[v].y, color);
		}
		DrawLineSDL(screen, projPos[0], projPos[1], color);
		DrawLineSDL(screen, projPos[1], projPos[2], color);
		DrawLineSDL(screen, projPos[0], projPos[2], color);*/


	}



	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

void VertexShader(const vec3& v, ivec2& p)
{

	vec3 v_trans = v - cameraPos;
	float f = SCREEN_WIDTH / 1;
	p.x = f*v_trans.x / v_trans.z + SCREEN_WIDTH / 2;
	p.y = f*v_trans.y / v_trans.z + SCREEN_HEIGHT / 2;

}

void Interpolate(ivec2 a, ivec2 b, vector<ivec2>& result)
{
	int N = result.size();
	vec2 step = vec2(b - a) / float(fmax(N - 1, 1));
	vec2 current(a);
	vec2 current_round=current;
	for (int i = 0; i<N; ++i)
	{
		result[i] = current_round;
		current += step;
		current_round.x = round(current.x);
		current_round.y = round(current.y);

	}
}

void DrawLineSDL(SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color)
{
	ivec2 delta = glm::abs(a - b);
	int pixels = glm::max(delta.x, delta.y) + 1;
	vector<ivec2> line(pixels);
	Interpolate(a, b, line);
	for (int j = 0; j < line.size(); ++j)
	{
		PutPixelSDL(screen, line[j].x, line[j].y, color);
	}
}

void DrawPolygonEdges(const vector<vec3>& vertices)
{
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position: 
	vector<ivec2> projectedVertices(V);
	for (int i = 0; i<V; ++i)
	{
		VertexShader(vertices[i], projectedVertices[i]);
	}


	// Loop over all vertices and draw the edge from it to the next vertex: 
	for (int i = 0; i<V; ++i)
	{
		int j = (i + 1) % V; // The next vertex 
		vec3 color(1, 1, 1);
		DrawLineSDL(screen, projectedVertices[i], projectedVertices[j],
			color);
	}
}


void ComputePolygonRows(const vector<ivec2>& vertexPixels, vector<ivec2>& leftPixels,
	vector<ivec2>& rightPixels)
{
	// 1. Find max and min y-value of the polygon 
	// and compute the number of rows it occupies. 
	int y_max = -numeric_limits<int>::max();
	int y_min = numeric_limits<int>::max();
	for (int i = 0; i < vertexPixels.size(); ++i)
	{
		if (vertexPixels[i].y > y_max) y_max = vertexPixels[i].y;
		if (vertexPixels[i].y < y_min) y_min = vertexPixels[i].y;
	}
	int rows = y_max - y_min + 1;

	// 2. Resize leftPixels and rightPixels 
	// so that they have an element for each row. 

	leftPixels.resize(rows);
	rightPixels.resize(rows);

	// 3. Initialize the x-coordinates in leftPixels 
	// to some really large value and the x-coordinates 
	// in rightPixels to some really small value. 

	for (int i = 0; i<rows; ++i)
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
	}

	// 4. Loop through all edges of the polygon and use 
	// linear interpolation to find the x-coordinate for 
	// each row it occupies. Update the corresponding 
	// values in rightPixels and leftPixels. 


	for (int i = 0; i < vertexPixels.size() - 1; ++i)
	{
		for (int j = i + 1; j < vertexPixels.size(); j++)
		{
			vector<ivec2> lines;
			ivec2 delta = glm::abs(vertexPixels[i] - vertexPixels[j]);
			int pixels = glm::max(delta.x, delta.y) + 1;
			lines.resize(pixels);
			Interpolate(vertexPixels[i], vertexPixels[j], lines);
			int seq;
			for (int i = 0; i < lines.size(); ++i)
			{
				seq = lines[i].y - y_min;
				if (leftPixels[seq].x > lines[i].x) 
				{ 
					leftPixels[seq] = lines[i]; 
				}
				if (rightPixels[seq].x < lines[i].x)
				{
					rightPixels[seq] = lines[i];
				}
			}
		}
	}

}

void DrawRows(const vector<ivec2>& leftPixels,
	const vector<ivec2>& rightPixels)
{
	int N = leftPixels.size();
	for (int y = 0; y < N; ++y)
	{
		for (int x = leftPixels[y].x; x <= rightPixels[y].x; x++)
		{
			PutPixelSDL(screen, x, leftPixels[y].y, currentColor);
		}
	}
}

void DrawPolygon(const vector<Vertex>& vertices)
{
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);//vector<ivec2> vertexPixels(V);
	for (int i = 0; i < V; ++i)
	{
		VertexShader(vertices[i], vertexPixels[i]);
	}
	/*vector<ivec2> leftPixels;
	vector<ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawRows(leftPixels, rightPixels);*/

	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawPolygonRows(leftPixels, rightPixels);
}

//=============================================================
//===========================================================

void Interpolate(Pixel a, Pixel b, vector<Pixel>& result)
{
	int N = result.size();
	float step_x,step_y,step_zinv;
	vec3 step_ilu, ilu;
	step_ilu = (b.illumination - a.illumination) / float(fmax(N - 1, 1));
	step_x= (b.x - a.x) / float(fmax(N - 1, 1));
	step_y = (b.y - a.y) / float(fmax(N - 1, 1));
	step_zinv = (b.zinv - a.zinv) / float(fmax(N - 1, 1));
	vec2 current;
	current.x = a.x;
	current.y = a.y;
	Pixel current_round=a;
	for (int i = 0; i<N; ++i)
	{
		result[i] = current_round;
		current.x += step_x;
		current.y += step_y;
		current_round.zinv += step_zinv;
		current_round.illumination += step_ilu;
		current_round.x = round(current.x);
		current_round.y = round(current.y);
	}
}

void ComputePolygonRows(const vector<Pixel>& vertexPixels,
	vector<Pixel>& leftPixels,
	vector<Pixel>& rightPixels
	)
{
	// 1. Find max and min y-value of the polygon 
	// and compute the number of rows it occupies. 
	int y_max = -numeric_limits<int>::max();
	int y_min = numeric_limits<int>::max();
	for (int i = 0; i < vertexPixels.size(); ++i)
	{
		if (vertexPixels[i].y > y_max) y_max = vertexPixels[i].y;
		if (vertexPixels[i].y < y_min) y_min = vertexPixels[i].y;
	}
	int rows = y_max - y_min + 1;

	// 2. Resize leftPixels and rightPixels 
	// so that they have an element for each row. 

	leftPixels.resize(rows);
	rightPixels.resize(rows);

	// 3. Initialize the x-coordinates in leftPixels 
	// to some really large value and the x-coordinates 
	// in rightPixels to some really small value. 

	for (int i = 0; i<rows; ++i)
	{
		leftPixels[i].x = +numeric_limits<int>::max();
		rightPixels[i].x = -numeric_limits<int>::max();
	}

	// 4. Loop through all edges of the polygon and use 
	// linear interpolation to find the x-coordinate for 
	// each row it occupies. Update the corresponding 
	// values in rightPixels and leftPixels. 


	for (int i = 0; i < vertexPixels.size() - 1; ++i)
	{
		for (int j = i + 1; j < vertexPixels.size(); j++)
		{
			vector<Pixel> lines;
			Pixel delta;
			delta.x= abs(vertexPixels[i].x - vertexPixels[j].x);
			delta.y = abs(vertexPixels[i].y - vertexPixels[j].y);
			int pixels = glm::max(delta.x, delta.y) + 1;
			lines.resize(pixels);
			Interpolate(vertexPixels[i], vertexPixels[j], lines);
			int seq;
			for (int i = 0; i < lines.size(); ++i)
			{
				seq = lines[i].y - y_min;
				if (leftPixels[seq].x > lines[i].x)
				{
					leftPixels[seq]= lines[i];
				}
				if (rightPixels[seq].x < lines[i].x)
				{
					rightPixels[seq] = lines[i];
				}
			}
		}
	}
}

void DrawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels)
{
	
	int N = leftPixels.size();
	for (int y = 0; y < N; ++y)
	{
		int S = -leftPixels[y].x + rightPixels[y].x + 1;
		vector<Pixel> result(S);
		Interpolate(leftPixels[y], rightPixels[y], result);
		//for (int x = leftPixels[y].x; x <= rightPixels[y].x; ++x)
		//{
		//	PutPixelSDL(screen, x, leftPixels[y].y, currentColor);
		//}

		for (int i = 0; i <S; i++)
		{				
			PixelShader(result[i]);
		}
	}
}
void VertexShader(const Vertex& v, Pixel& p)
{
	vec3 v_trans = v.position - cameraPos;
	float f = SCREEN_WIDTH / 1;
	p.x = f*v_trans.x / v_trans.z + SCREEN_WIDTH / 2;
	p.y = f*v_trans.y / v_trans.z + SCREEN_HEIGHT / 2;
	p.zinv = 1.0 / v_trans.z;
	vec3 r = lightPos - v.position;
	vec3 D = lightPower*float(fmax(glm::dot(glm::normalize(r), v.normal), 0) / 4 * glm::dot(r, r)*3.1415926f);
	p.illumination = v.reflectance*(indirectLightPowerPerArea+D);
}


void PixelShader(const Pixel& p)
{
	int x = p.x;
	int y = p.y;
	if (p.zinv > depthBuffer[y][x])
	{
		depthBuffer[y][x] = p.zinv;
		PutPixelSDL(screen, x, y, p.illumination);
	}
}