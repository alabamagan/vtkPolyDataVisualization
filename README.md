*Usage* 

./vtkVisualization -i <vtkfile> [-c <camera>|-a <azimuth>|-e <elevation>|-s 400x400|-o <./tmp.jpg>]

*Options*:
	-h,--help		Show this help message
	-i,--input		Directory to input vtk data file.
	-a,--azimuth		Azimuthal rotation applied to the model along the axis that run in the direction of view up angle through the focal point
	-e,--elevation		Elevation rotation applied to the model.
	-c,--camera		Points to the file that stores previous camera information
	-o,--output		Output image directory, default to [./tmp.jpg]
