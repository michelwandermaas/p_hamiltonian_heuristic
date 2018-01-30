#include "Graph.h"

Graph::Graph()
{
	M = NULL;
	X = NULL;
	Y = NULL;
}

Graph::Graph(char *filename)
{
	Read(filename);
}

Graph::~Graph()
{
	if(M)
	{
		for(int i = 0; i < n; i++)
		{
			delete [] M[i];
		}
		delete [] M;
	}
	if(X)
	{
		delete X;
	}
	if(Y)
	{
		delete Y;
	}
}

void Graph::Read(char *filename)
{
	ifstream file;
	file.open(filename);

	if(!file.is_open())
	{
		printf("Couldnt open the file %s\n", filename);
		exit(0);
	}

	string name, type, EdgeWeightType, EdgeWeightFormat, EdgeDataType, NodeCoordType, DisplayDataType, line;

	X = NULL;
	Y = NULL;

	while(!file.eof())
	{
		string keyword;
		stringstream stream;
		
		getline(file, line);
		stream << line;		

		getline(stream, keyword, ':');
		size_t pos = keyword.find(" ");
		while(pos != string::npos)
		{
			keyword.erase(pos);
			pos = keyword.find(" ");
		}

		if(keyword == "NAME")
			stream >> name;
		else if(keyword == "TYPE")
			stream >> type;
		else if(keyword == "DIMENSION")
		{
			stream >> n;

			M = new int*[n];
			for(int i = 0; i < n; i++)
			{
				M[i] = new int[n];
				for(int j= 0; j < n; j++)
				{
					M[i][j] = -1;
				}
			}
		}
		else if(keyword == "EDGE_WEIGHT_TYPE")
			stream >> EdgeWeightType;
		else if(keyword == "EDGE_WEIGHT_FORMAT")
			stream >> EdgeWeightFormat;
		else if(keyword == "EDGE_DATA_TYPE")
			stream >> EdgeDataType;
		else if(keyword == "NODE_COORD_TYPE")
			stream >> NodeCoordType;
		else if(keyword == "DISPLAY_DATA_TYPE")
			stream >> DisplayDataType;
		else if(keyword == "NODE_COORD_SECTION")
		{
			if(NodeCoordType == "TWOD_COORDS" ||
				EdgeWeightType == "EUC_2D" ||
				EdgeWeightType == "MAX_2D" ||
				EdgeWeightType == "CEIL_2D" ||
				EdgeWeightType == "ATT" ||
				EdgeWeightType == "GEO")
			{
				X = new double[n];
				Y = new double[n];
				int id;
				for(int i = 0; i < n; i++)
				{
					file >> id >> X[i] >> Y[i];
					if(EdgeWeightType == "GEO")
					{
						double PI = 3.141592;
						double deg = round(X[i]);
						double min = X[i] - deg;
						X[i] = PI * ( deg + 5.0 * min / 3.0 ) / 180.0;
						deg = round(Y[i]);
						min = Y[i] - deg;
						Y[i] = PI * ( deg + 5.0 * min / 3.0 ) / 180.0;
					}
				}
				for(int i = 0; i < n; i++)
				{
					for(int j = 0; j < n; j++)
					{
						if(i==j) continue;
				
						if(EdgeWeightType == "EUC_2D")
						{
							double xd = X[i] - X[j];
							double yd = Y[i] - Y[j];
							M[i][j] = M[j][i] = round( sqrt( xd*xd + yd*yd ) );
						}
						else if(EdgeWeightType == "MAX_2D")
						{
							double xd = fabs( X[i] - X[j] );
							double yd = fabs( Y[i] - Y[j] );
							M[i][j] = M[j][i] = max( round(xd), round(yd) );
						}
						else if(EdgeWeightType == "GEO")
						{
							double RRR = 6378.388;
							double q1 = cos( Y[i] - Y[j] );
							double q2 = cos( X[i] - X[j] );
							double q3 = cos( X[i] + X[j] );
							M[i][j] = M[j][i] = (int)( RRR * acos( 0.5 * ( ( 1.0 + q1 ) * q2 - (1.0 - q1) * q3 ) ) + 1.0 );
						}
						else if(EdgeWeightType == "ATT")
						{
							double xd = X[i] - X[j];
							double yd = Y[i] - Y[j];
							double rij = sqrt( ( xd * xd + yd * yd ) / 10.0);
							double tij = round( rij );
							if(tij < rij ) M[i][j] = M[j][i] = tij + 1.0;
							else M[i][j] = M[j][i] = tij;
						}
						else if(EdgeWeightType == "CEIL_2D")
						{
							double xd = X[i] - X[j];
							double yd = Y[i] - Y[j];
							M[i][j] = M[j][i] = (int)( sqrt( xd*xd + yd*yd ) + 1.0);
						}
					}
				}
			}
			else
			{
				printf("data format not supported\n");
				exit(0);
			}
		}
		else if(keyword == "EDGE_WEIGHT_SECTION")
		{
			if(EdgeWeightFormat == "FULL_MATRIX")
			{
				for(int i = 0; i < n; i++)
				{
					for(int j = 0; j < n; j++)
					{
						file >> M[i][j];
					}
				}
			}
			else if(EdgeWeightFormat == "UPPER_ROW")
			{
				for(int i = 0; i < n; i++)
				{
					for(int j = i +1; j < n; j++)
					{
						file >> M[i][j];
						M[j][i] = M[i][j];
					}
				}
			}
			else if(EdgeWeightFormat == "LOWER_ROW")
			{
				for(int i = 0; i < n; i++)
				{
					for(int j = 0; j < i; j++)
					{
						file >> M[i][j];
						M[j][i] = M[i][j];
					}
				}
			}
			else if(EdgeWeightFormat == "UPPER_DIAG_ROW")
			{
				for(int i = 0; i < n; i++)
				{
					int t;
					file >> t;
					for(int j = i+1; j < n; j++)
					{
						file >> M[i][j];
						M[j][i] = M[i][j];
					}
				}
			}
			else if(EdgeWeightFormat == "LOWER_DIAG_ROW")
			{
				for(int i = 0; i < n; i++)
				{
					for(int j = 0; j < i; j++)
					{
						file >> M[i][j];
						M[j][i] = M[i][j];
					}
					int t;
					file >> t;
				}
			}
			else if(EdgeWeightFormat == "UPPER_COL")
			{
				printf("data format not supported\n");
				exit(0);
			}
			else if(EdgeWeightFormat == "LOWER_COL")
			{
				printf("data format not supported\n");
				exit(0);
			}
			else if(EdgeWeightFormat == "UPPER_DIAG_COL")
			{
				printf("data format not supported\n");
				exit(0);
			}
			else if(EdgeWeightFormat == "LOWER_DIAG_COL")
			{
				printf("data format not supported\n");
				exit(0);
			}
		}
		else if(keyword == "DISPLAY_DATA_SECTION")
		{
			if(DisplayDataType != "TWOD_COORDS") continue;

			if(!X || !Y)
			{
				X = new double[n];
				Y = new double [n];
			}
			for(int i = 0; i < n; i++)
			{
				int id;
				file >> id >> X[i] >> Y[i];
			}
		}
		else
		{
			getline(stream, keyword);
			//file.getline((char*)&ch, 100);
			//cout << keyword << endl;
		}
	}
	
	file.close();
}

