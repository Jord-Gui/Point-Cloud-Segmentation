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