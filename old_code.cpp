/*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
OLD BOUNDING BOX FINDER - VISUALISATION CAUSES ERRORS 13/10/2019*/
//#include <iostream>
//#include <vector>
//#include <pcl/point_types.h>
//#include <pcl/io/pcd_io.h>
//#include <pcl/search/search.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/filters/passthrough.h>
//#include <pcl/segmentation/region_growing.h>
//
//#include <pcl/common/common.h>
//#include <pcl/common/pca.h>
//#include <pcl/common/transforms.h>
//
//int
//main(int argc, char** argv)
//{
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
//	if (pcl::io::loadPCDFile <pcl::PointXYZ>("region_growing_tutorial.pcd", *cloud) == -1)
//	{
//		std::cout << "Cloud reading failed." << std::endl;
//		return (-1);
//	}
//
//	pcl::search::Search<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
//	pcl::PointCloud <pcl::Normal>::Ptr normals(new pcl::PointCloud <pcl::Normal>);
//	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normal_estimator;
//	normal_estimator.setSearchMethod(tree);
//	normal_estimator.setInputCloud(cloud);
//	normal_estimator.setKSearch(50);
//	normal_estimator.compute(*normals);
//
//	pcl::IndicesPtr indices(new std::vector <int>);
//	pcl::PassThrough<pcl::PointXYZ> pass;
//	pass.setInputCloud(cloud);
//	pass.setFilterFieldName("z");
//	pass.setFilterLimits(0.0, 1.0);
//	pass.filter(*indices);
//
//	pcl::RegionGrowing<pcl::PointXYZ, pcl::Normal> reg;
//	reg.setMinClusterSize(50);
//	reg.setMaxClusterSize(1000000);
//	reg.setSearchMethod(tree);
//	reg.setNumberOfNeighbours(30);
//	reg.setInputCloud(cloud);
//	//reg.setIndices (indices);
//	reg.setInputNormals(normals);
//	reg.setSmoothnessThreshold(3.0 / 180.0 * M_PI);
//	reg.setCurvatureThreshold(1.0);
//
//	std::vector <pcl::PointIndices> clusters;
//	reg.extract(clusters);
//
//	pcl::PointCloud <pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud();
//
//	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloudPCAprojection(new pcl::PointCloud<pcl::PointXYZRGB>);
//	pcl::PCA<pcl::PointXYZRGB> pca;
//	pca.setInputCloud(colored_cloud);
//	pca.project(*colored_cloud, *cloudPCAprojection);
//	std::cerr << std::endl << "EigenVectors: " << pca.getEigenVectors() << std::endl;
//	std::cerr << std::endl << "EigenValues: " << pca.getEigenValues() << std::endl;
//
//	// Transform the original cloud to the origin where the principal components correspond to the axes.
//	Eigen::Matrix4f projectionTransform(Eigen::Matrix4f::Identity());
//	Eigen::Vector4f pcaCentroid;
//	Eigen::Matrix3f covariance;
//	computeCovarianceMatrixNormalized(*colored_cloud, pcaCentroid, covariance);
//	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
//	Eigen::Matrix3f eigenVectorsPCA = eigen_solver.eigenvectors();
//	projectionTransform.block<3, 3>(0, 0) = eigenVectorsPCA.transpose();
//	projectionTransform.block<3, 1>(0, 3) = -1.f * (projectionTransform.block<3, 3>(0, 0) * pcaCentroid.head<3>());
//	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloudPointsProjected(new pcl::PointCloud<pcl::PointXYZRGB>);
//	pcl::transformPointCloud(*colored_cloud, *cloudPointsProjected, projectionTransform);
//	// Get the minimum and maximum points of the transformed cloud.
//	pcl::PointXYZRGB minPoint, maxPoint;
//	//need #include <pcl/common/common.h>
//	pcl::getMinMax3D(*cloudPointsProjected, minPoint, maxPoint);
//	const Eigen::Vector3f meanDiagonal = 0.5f*(maxPoint.getVector3fMap() + minPoint.getVector3fMap());
//
//	// Final transform
//	const Eigen::Quaternionf bboxQuaternion(eigenVectorsPCA); //Quaternions are a way to do rotations
//	const Eigen::Vector3f bboxTransform = eigenVectorsPCA * meanDiagonal + pcaCentroid.head<3>();
//
//	std::cout << "Before Visualisation" << std::endl;
//
//	// This viewer has 4 windows, but is only showing images in one of them as written here.
//	pcl::visualization::PCLVisualizer *visu;
//	visu = new pcl::visualization::PCLVisualizer(argc, argv, "PlyViewer");
//	int mesh_vp_1, mesh_vp_2, mesh_vp_3, mesh_vp_4;
//	visu->createViewPort(0.0, 0.5, 0.5, 1.0, mesh_vp_1);
//	visu->createViewPort(0.5, 0.5, 1.0, 1.0, mesh_vp_2);
//	visu->createViewPort(0.0, 0, 0.5, 0.5, mesh_vp_3);
//	visu->createViewPort(0.5, 0, 1.0, 0.5, mesh_vp_4);
//	visu->addPointCloud(colored_cloud, "bboxedCloud", mesh_vp_3);
//	visu->addCube(bboxTransform, bboxQuaternion, maxPoint.x - minPoint.x, maxPoint.y - minPoint.y, maxPoint.z - minPoint.z, "bbox", mesh_vp_3);
//
//	return (0);
//}


/*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
SEGMENTATION ALGORITHM THAT OUTPUTS EACH CLUSTER INTO THEIR OWN SEPERATE .PCD FILE 15/10/2019*/
//#include <iostream>
//#include <vector>
//#include <pcl/point_types.h>
//#include <pcl/io/pcd_io.h>
//#include <pcl/search/search.h>
//#include <pcl/search/kdtree.h>
//#include <pcl/features/normal_3d.h>
//#include <pcl/visualization/cloud_viewer.h>
//#include <pcl/filters/passthrough.h>
//#include <pcl/segmentation/region_growing.h>
//
//int
//main(int argc, char** argv)
//{
//	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
//	if (pcl::io::loadPCDFile <pcl::PointXYZ>("region_growing_tutorial.pcd", *cloud) == -1)
//	{
//		std::cout << "Cloud reading failed." << std::endl;
//		return (-1);
//	}
//
//	pcl::search::Search<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
//	pcl::PointCloud <pcl::Normal>::Ptr normals(new pcl::PointCloud <pcl::Normal>);
//	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normal_estimator;
//	normal_estimator.setSearchMethod(tree);
//	normal_estimator.setInputCloud(cloud);
//	normal_estimator.setKSearch(50);
//	normal_estimator.compute(*normals);
//
//	pcl::IndicesPtr indices(new std::vector <int>);
//	pcl::PassThrough<pcl::PointXYZ> pass;
//	pass.setInputCloud(cloud);
//	pass.setFilterFieldName("z");
//	pass.setFilterLimits(0.0, 1.0);
//	pass.filter(*indices);
//
//	pcl::RegionGrowing<pcl::PointXYZ, pcl::Normal> reg;
//	reg.setMinClusterSize(50);
//	reg.setMaxClusterSize(1000000);
//	reg.setSearchMethod(tree);
//	reg.setNumberOfNeighbours(30);
//	reg.setInputCloud(cloud);
//	//reg.setIndices (indices);
//	reg.setInputNormals(normals);
//	reg.setSmoothnessThreshold(3.0 / 180.0 * M_PI);
//	reg.setCurvatureThreshold(1.0);
//
//	std::vector <pcl::PointIndices> clusters;
//	reg.extract(clusters);
//
//	// my code
//	std::cout << "Segmentation done!" << std::endl;
//	
//	pcl::PointCloud <pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud();
//	int counter = 0;
//	while (counter < clusters.size()) {
//		pcl::PointCloud <pcl::PointXYZRGB>::Ptr extracted_cloud(new pcl::PointCloud <pcl::PointXYZRGB>);
//		pcl::copyPointCloud(*colored_cloud, clusters[counter], *extracted_cloud);
//		std::string fileName = "cluster_";
//		fileName += std::to_string(counter) + ".pcd";
//		pcl::io::savePCDFileASCII(fileName, *extracted_cloud);
//		std::cout << "Finished!" << std::endl;
//		counter++;
//	}
//
//	return (0);
//}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
FIND BOUNDING BOX OF A SINGLE POINT CLOUD AND VISUALISATION WORKS 15/10/2019*/
//#include <vtkAutoInit.h>         
//VTK_MODULE_INIT(vtkRenderingOpenGL);
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
//
//#include <iostream>
//#include <string>
//#include <pcl/io/pcd_io.h>
//#include <pcl/point_cloud.h>
//#include <pcl/point_types.h>
//#include <Eigen/Core>
//#include <pcl/common/transforms.h>
//#include <pcl/common/common.h>
//#include <pcl/visualization/pcl_visualizer.h>
//
//using namespace std;
//typedef pcl::PointXYZ PointType;
//
//int main(int argc, char **argv)
//{
//	pcl::PointCloud<PointType>::Ptr cloud(new pcl::PointCloud<PointType>());
//
//	pcl::io::loadPCDFile("cluster_71.pcd", *cloud);
//
//	/*
//	The pca principal component analysis method is used to obtain the three main directions of the point cloud,
//	obtain the centroid, calculate the covariance, obtain the covariance matrix,
//	and obtain the eigenvalue and feature vector of the covariance matrix.
//	the feature vector is the main direction.
//	*/
//
//	Eigen::Vector4f pcaCentroid;
//	pcl::compute3DCentroid(*cloud, pcaCentroid);
//	Eigen::Matrix3f covariance;
//	pcl::computeCovarianceMatrixNormalized(*cloud, pcaCentroid, covariance);
//	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
//	Eigen::Matrix3f eigenVectorsPCA = eigen_solver.eigenvectors();
//	Eigen::Vector3f eigenValuesPCA = eigen_solver.eigenvalues();
//	eigenVectorsPCA.col(2) = eigenVectorsPCA.col(0).cross(eigenVectorsPCA.col(1)); // Correct vertical between main directions
//	eigenVectorsPCA.col(0) = eigenVectorsPCA.col(1).cross(eigenVectorsPCA.col(2));
//	eigenVectorsPCA.col(1) = eigenVectorsPCA.col(2).cross(eigenVectorsPCA.col(0));
//
//	std::cout << "Eigenvalues (3x1):\n" << eigenValuesPCA << std::endl;
//	std::cout << "Eigenvectors (3x3):\n" << eigenVectorsPCA << std::endl;
//	std::cout << "Centroid point (4x1):\n" << pcaCentroid << std::endl;
//
//	/*
//	Using the main direction and centroid obtained above, the input point cloud is converted to the origin,
//	and the main direction. the coordinate system direction is returned,
//	and a bounding box of the point cloud transformed to the origin is established.
//	*/
//
//	Eigen::Matrix4f tm = Eigen::Matrix4f::Identity();
//	Eigen::Matrix4f tm_inv = Eigen::Matrix4f::Identity();
//	tm.block<3, 3>(0, 0) = eigenVectorsPCA.transpose();   //R.
//	tm.block<3, 1>(0, 3) = -1.0f * (eigenVectorsPCA.transpose()) *(pcaCentroid.head<3>());//  -R*t
//	tm_inv = tm.inverse();
//
//	std::cout << "Transformation matrix (4x4):\n" << tm << std::endl;
//	std::cout << "Transformation matrix (4x4):\n" << tm_inv << std::endl;
//
//	pcl::PointCloud<PointType>::Ptr transformedCloud(new pcl::PointCloud<PointType>);
//	pcl::transformPointCloud(*cloud, *transformedCloud, tm);
//
//	PointType min_p1, max_p1;
//	Eigen::Vector3f c1, c;
//	pcl::getMinMax3D(*transformedCloud, min_p1, max_p1);
//	c1 = 0.5f*(min_p1.getVector3fMap() + max_p1.getVector3fMap());
//
//	std::cout << "Heart c1(3x1):\n" << c1 << std::endl;
//
//	Eigen::Affine3f tm_inv_aff(tm_inv);
//	pcl::transformPoint(c1, c, tm_inv_aff);
//
//	Eigen::Vector3f whd, whd1;
//	whd1 = max_p1.getVector3fMap() - min_p1.getVector3fMap();
//	whd = whd1;
//	float sc1 = (whd1(0) + whd1(1) + whd1(2)) / 3;  // Point cloud average scale, used to set the main direction arrow size
//
//	std::cout << "width1=" << whd1(0) << endl;
//	std::cout << "heght1=" << whd1(1) << endl;
//	std::cout << "depth1=" << whd1(2) << endl;
//	std::cout << "scale1=" << sc1 << endl;
//
//	const Eigen::Quaternionf bboxQ1(Eigen::Quaternionf::Identity());
//	const Eigen::Vector3f    bboxT1(c1);
//
//	const Eigen::Quaternionf bboxQ(tm_inv.block<3, 3>(0, 0));
//	const Eigen::Vector3f    bboxT(c);
//
//	/*
//	Set the main direction and bounding box to the input point cloud,
//	and implement the inverse transformation by inputting the point cloud to the origin point cloud transform.
//	*/
//
//	// The main direction of the point cloud that is transformed to the origin
//	PointType op;
//	op.x = 0.0;
//	op.y = 0.0;
//	op.z = 0.0;
//	Eigen::Vector3f px, py, pz;
//	Eigen::Affine3f tm_aff(tm);
//	pcl::transformVector(eigenVectorsPCA.col(0), px, tm_aff);
//	pcl::transformVector(eigenVectorsPCA.col(1), py, tm_aff);
//	pcl::transformVector(eigenVectorsPCA.col(2), pz, tm_aff);
//	PointType pcaX;
//	pcaX.x = sc1 * px(0);
//	pcaX.y = sc1 * px(1);
//	pcaX.z = sc1 * px(2);
//	PointType pcaY;
//	pcaY.x = sc1 * py(0);
//	pcaY.y = sc1 * py(1);
//	pcaY.z = sc1 * py(2);
//	PointType pcaZ;
//	pcaZ.x = sc1 * pz(0);
//	pcaZ.y = sc1 * pz(1);
//	pcaZ.z = sc1 * pz(2);
//
//	// The main direction of the initial point cloud
//	PointType cp;
//	cp.x = pcaCentroid(0);
//	cp.y = pcaCentroid(1);
//	cp.z = pcaCentroid(2);
//	PointType pcX;
//	pcX.x = sc1 * eigenVectorsPCA(0, 0) + cp.x;
//	pcX.y = sc1 * eigenVectorsPCA(1, 0) + cp.y;
//	pcX.z = sc1 * eigenVectorsPCA(2, 0) + cp.z;
//	PointType pcY;
//	pcY.x = sc1 * eigenVectorsPCA(0, 1) + cp.x;
//	pcY.y = sc1 * eigenVectorsPCA(1, 1) + cp.y;
//	pcY.z = sc1 * eigenVectorsPCA(2, 1) + cp.z;
//	PointType pcZ;
//	pcZ.x = sc1 * eigenVectorsPCA(0, 2) + cp.x;
//	pcZ.y = sc1 * eigenVectorsPCA(1, 2) + cp.y;
//	pcZ.z = sc1 * eigenVectorsPCA(2, 2) + cp.z;
//
//	// Visualization
//	pcl::visualization::PCLVisualizer viewer;
//
//	pcl::visualization::PointCloudColorHandlerCustom<PointType> tc_handler(transformedCloud, 0, 255, 0); // Point cloud related to the origin
//	viewer.addPointCloud(transformedCloud, tc_handler, "transformCloud");
//	viewer.addCube(bboxT1, bboxQ1, whd1(0), whd1(1), whd1(2), "bbox1");
//	viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME, "bbox1");
//	viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "bbox1");
//
//	pcl::visualization::PointCloudColorHandlerCustom<PointType> color_handler(cloud, 255, 0, 0);  // Input initial point cloud related
//	viewer.addPointCloud(cloud, color_handler, "cloud");
//	viewer.addCube(bboxT, bboxQ, whd(0), whd(1), whd(2), "bbox");
//	viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME, "bbox");
//	viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "bbox");
//
//	viewer.addCoordinateSystem(0.5f*sc1); // Adds a xyz axis at (0, 0, 0) of viewer
//	viewer.setBackgroundColor(1.0, 1.0, 1.0);
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//	}
//
//	return 0;
//}

/*------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
BOUNDING BOX FOR MULTIPLE FILES IN VIEWER BUT CANNOT SAVE THE VIEWER 18/10/2019*/
//#include <vtkAutoInit.h>         
//VTK_MODULE_INIT(vtkRenderingOpenGL);
//VTK_MODULE_INIT(vtkInteractionStyle);
//VTK_MODULE_INIT(vtkRenderingFreeType);
//
//#include <iostream>
//#include <string>
//#include <pcl/io/pcd_io.h>
//#include <pcl/point_cloud.h>
//#include <pcl/point_types.h>
//#include <Eigen/Core>
//#include <pcl/common/transforms.h>
//#include <pcl/common/common.h>
//#include <pcl/visualization/pcl_visualizer.h>
//
//using namespace std;
//typedef pcl::PointXYZ PointType;
//
//int main(int argc, char **argv)
//{
//	int counter = 0;
//	int files = 133;
//	pcl::visualization::PCLVisualizer viewer;
//
//	while (counter < files) {
//		pcl::PointCloud<PointType>::Ptr cloud(new pcl::PointCloud<PointType>());
//
//		std::string fileName = "cluster_";
//		fileName += std::to_string(counter) + ".pcd";
//		pcl::io::loadPCDFile(fileName, *cloud);
//
//		/*
//		The pca principal component analysis method is used to obtain the three main directions of the point cloud,
//		obtain the centroid, calculate the covariance, obtain the covariance matrix,
//		and obtain the eigenvalue and feature vector of the covariance matrix.
//		the feature vector is the main direction.
//		*/
//
//		Eigen::Vector4f pcaCentroid;
//		pcl::compute3DCentroid(*cloud, pcaCentroid);
//		Eigen::Matrix3f covariance;
//		pcl::computeCovarianceMatrixNormalized(*cloud, pcaCentroid, covariance);
//		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
//		Eigen::Matrix3f eigenVectorsPCA = eigen_solver.eigenvectors();
//		Eigen::Vector3f eigenValuesPCA = eigen_solver.eigenvalues();
//		eigenVectorsPCA.col(2) = eigenVectorsPCA.col(0).cross(eigenVectorsPCA.col(1)); // Correct vertical between main directions
//		eigenVectorsPCA.col(0) = eigenVectorsPCA.col(1).cross(eigenVectorsPCA.col(2));
//		eigenVectorsPCA.col(1) = eigenVectorsPCA.col(2).cross(eigenVectorsPCA.col(0));
//
//		std::cout << "Eigenvalues (3x1):\n" << eigenValuesPCA << std::endl;
//		std::cout << "Eigenvectors (3x3):\n" << eigenVectorsPCA << std::endl;
//		std::cout << "Centroid point (4x1):\n" << pcaCentroid << std::endl;
//
//		/*
//		Using the main direction and centroid obtained above, the input point cloud is converted to the origin,
//		and the main direction. the coordinate system direction is returned,
//		and a bounding box of the point cloud transformed to the origin is established.
//		*/
//
//		Eigen::Matrix4f tm = Eigen::Matrix4f::Identity();
//		Eigen::Matrix4f tm_inv = Eigen::Matrix4f::Identity();
//		tm.block<3, 3>(0, 0) = eigenVectorsPCA.transpose();   //R.
//		tm.block<3, 1>(0, 3) = -1.0f * (eigenVectorsPCA.transpose()) *(pcaCentroid.head<3>());//  -R*t
//		tm_inv = tm.inverse();
//
//		std::cout << "Transformation matrix (4x4):\n" << tm << std::endl;
//		std::cout << "Transformation matrix (4x4):\n" << tm_inv << std::endl;
//
//		pcl::PointCloud<PointType>::Ptr transformedCloud(new pcl::PointCloud<PointType>);
//		pcl::transformPointCloud(*cloud, *transformedCloud, tm);
//
//		PointType min_p1, max_p1;
//		Eigen::Vector3f c1, c;
//		pcl::getMinMax3D(*transformedCloud, min_p1, max_p1);
//		c1 = 0.5f*(min_p1.getVector3fMap() + max_p1.getVector3fMap());
//
//		std::cout << "Heart c1(3x1):\n" << c1 << std::endl;
//
//		Eigen::Affine3f tm_inv_aff(tm_inv);
//		pcl::transformPoint(c1, c, tm_inv_aff);
//
//		Eigen::Vector3f whd, whd1;
//		whd1 = max_p1.getVector3fMap() - min_p1.getVector3fMap();
//		whd = whd1;
//		float sc1 = (whd1(0) + whd1(1) + whd1(2)) / 3;  // Point cloud average scale, used to set the main direction arrow size
//
//		std::cout << "width1=" << whd1(0) << endl;
//		std::cout << "heght1=" << whd1(1) << endl;
//		std::cout << "depth1=" << whd1(2) << endl;
//		std::cout << "scale1=" << sc1 << endl;
//
//		const Eigen::Quaternionf bboxQ1(Eigen::Quaternionf::Identity());
//		const Eigen::Vector3f    bboxT1(c1);
//
//		const Eigen::Quaternionf bboxQ(tm_inv.block<3, 3>(0, 0));
//		const Eigen::Vector3f    bboxT(c);
//
//
//		/*
//		Set the main direction and bounding box to the input point cloud,
//		and implement the inverse transformation by inputting the point cloud to the origin point cloud transform.
//		*/
//
//		// The main direction of the point cloud that is transformed to the origin
//		PointType op;
//		op.x = 0.0;
//		op.y = 0.0;
//		op.z = 0.0;
//		Eigen::Vector3f px, py, pz;
//		Eigen::Affine3f tm_aff(tm);
//		pcl::transformVector(eigenVectorsPCA.col(0), px, tm_aff);
//		pcl::transformVector(eigenVectorsPCA.col(1), py, tm_aff);
//		pcl::transformVector(eigenVectorsPCA.col(2), pz, tm_aff);
//		PointType pcaX;
//		pcaX.x = sc1 * px(0);
//		pcaX.y = sc1 * px(1);
//		pcaX.z = sc1 * px(2);
//		PointType pcaY;
//		pcaY.x = sc1 * py(0);
//		pcaY.y = sc1 * py(1);
//		pcaY.z = sc1 * py(2);
//		PointType pcaZ;
//		pcaZ.x = sc1 * pz(0);
//		pcaZ.y = sc1 * pz(1);
//		pcaZ.z = sc1 * pz(2);
//
//
//		// The main direction of the initial point cloud
//		PointType cp;
//		cp.x = pcaCentroid(0);
//		cp.y = pcaCentroid(1);
//		cp.z = pcaCentroid(2);
//		PointType pcX;
//		pcX.x = sc1 * eigenVectorsPCA(0, 0) + cp.x;
//		pcX.y = sc1 * eigenVectorsPCA(1, 0) + cp.y;
//		pcX.z = sc1 * eigenVectorsPCA(2, 0) + cp.z;
//		PointType pcY;
//		pcY.x = sc1 * eigenVectorsPCA(0, 1) + cp.x;
//		pcY.y = sc1 * eigenVectorsPCA(1, 1) + cp.y;
//		pcY.z = sc1 * eigenVectorsPCA(2, 1) + cp.z;
//		PointType pcZ;
//		pcZ.x = sc1 * eigenVectorsPCA(0, 2) + cp.x;
//		pcZ.y = sc1 * eigenVectorsPCA(1, 2) + cp.y;
//		pcZ.z = sc1 * eigenVectorsPCA(2, 2) + cp.z;
//
//
//		// Visualization
//		std::string boxId = "bbox";
//		boxId += std::to_string(counter);
//		viewer.addCube(bboxT, bboxQ, whd(0), whd(1), whd(2), boxId);
//		viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_REPRESENTATION, pcl::visualization::PCL_VISUALIZER_REPRESENTATION_WIREFRAME, boxId);
//		viewer.setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, boxId);
//
//		counter++;
//	}
//
//	viewer.setBackgroundColor(1.0, 1.0, 1.0);
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(100);
//	}
//
//	return 0;
//}
