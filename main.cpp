#include <string>
#include <fstream>
#include <filesystem>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/utilities/eigen_interop_helpers.h>
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage);

Eigen::Map<Eigen::Vector3d> GCVec3toEigenVec3(geometrycentral::Vector3& vector);
Eigen::Map<Eigen::VectorXd> FlattenMatrixN3(Eigen::Matrix<double, 3, Eigen::Dynamic>& vector);
size_t pointToIndex(const geometrycentral::Vector3& point,
                  geometrycentral::surface::VertexPositionGeometry& geometry,
                  geometrycentral::surface::SurfaceMesh& mesh);


float structSpringConst = 35.0f;
float bendingSpringConst = 3.5f;
float dampingConstant = 0.1;
float massValue = 0.1;
float timeStep = 1.0f/60;
Eigen::Vector3d gravity = { 0, -9.8 , 0};

int main(int argc, char *argv[])
{
    int numberOfFiles;
    if ( argc == 1 )
    {
        numberOfFiles = 200;
    }
    else
    {
        numberOfFiles = std::atoi(argv[1]);
    }

    std::unique_ptr<geometrycentral::surface::SurfaceMesh> me;
	std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> geo;
	std::tie(me, geo) = geometrycentral::surface::readSurfaceMesh("../assets/square.obj");

    geometrycentral::surface::VertexPositionGeometry& geometry = *geo;
    geometrycentral::surface::SurfaceMesh& mesh = *me;

    Eigen::Matrix<double, 3, Eigen::Dynamic> forces;
    Eigen::Matrix<double, 3, Eigen::Dynamic> invMass;
    Eigen::VectorXd mass;
    Eigen::Matrix<double, 3, Eigen::Dynamic> velocities;
    Eigen::Map<Eigen::VectorXd, 16> positions = FlattenedEigenMap<double, 3>(geometry.vertexPositions);
    geometrycentral::surface::EdgeData<double> structLengths;
	geometrycentral::surface::EdgeData<double> bendingLengths;


    mass.resize(geometry.vertexPositions.size());
    forces.resize(3, geometry.vertexPositions.size());
    invMass.resize(3, geometry.vertexPositions.size());
    velocities.resize(3, geometry.vertexPositions.size());
    forces.setZero();
    velocities.setZero();

    geometry.requireVertexIndices();

    for (geometrycentral::surface::Vertex v : mesh.vertices())
    {
        size_t i = geometry.vertexIndices[v];
        invMass.col(i) = Eigen::Vector3d::Constant(1.0 / massValue);
        mass(i) = massValue;
    }

    // Fixing two points
    invMass.col(pointToIndex({-1.0, 0.0, 1.0}, geometry, mesh)) = Eigen::Vector3d::Constant(0.0);
    //geometry.vertexPositions[pointToIndex({-1.0, 0.0, 1.0}, geometry, mesh)] = {-1.0, 0.0, 1.0};

    invMass.col(pointToIndex({1.0, 0.0, 1.0}, geometry, mesh)) = Eigen::Vector3d::Constant(0.0);
    //geometry.vertexPositions[pointToIndex({1.0, 0.0, 1.0}, geometry, mesh)] = {1.0, 0.0, 1.0};

    // Calculating initial data
    structLengths = geometrycentral::surface::EdgeData<double>(mesh);
    for (geometrycentral::surface::Edge e : mesh.edges())
    {
        structLengths[e] = (geometry.vertexPositions[e.firstVertex()] - geometry.vertexPositions[e.secondVertex()]).norm();
    }


    bendingLengths = geometrycentral::surface::EdgeData<double>(mesh);
    for (geometrycentral::surface::Edge e : mesh.edges())
    {
        geometrycentral::surface::Halfedge he = e.halfedge();
        geometrycentral::surface::Vertex vertex1 = he.next().next().vertex();
        geometrycentral::surface::Vertex vertex2 = he.twin().next().next().vertex();

        bendingLengths[e] = (geometry.vertexPositions[vertex1] - geometry.vertexPositions[vertex2]).norm();
    }

    // Running the simulation
    for (int frame = 0; frame < numberOfFiles; frame++)
    {
        // Progress Bar
        printProgress(float(frame + 1)/numberOfFiles);

        geometry.requireVertexIndices();
        forces.setZero();
        // Adding Structural Springs
        for (geometrycentral::surface::Edge e : mesh.edges())
        {
            geometrycentral::surface::Vertex vertex1 = e.firstVertex();
            geometrycentral::surface::Vertex vertex2 = e.secondVertex();

            const size_t index1 = geometry.vertexIndices[vertex1];
            const size_t index2 = geometry.vertexIndices[vertex2];

            const Eigen::Vector3d diffVector = GCVec3toEigenVec3(geometry.vertexPositions[index2]) - GCVec3toEigenVec3(geometry.vertexPositions[index1]);
            const Eigen::Vector3d direction = diffVector.normalized();
            const double forceMagnitude = (diffVector.norm() - structLengths[e]) * structSpringConst;
            const Eigen::Vector3d force = forceMagnitude * direction;

            const Eigen::Vector3d& velocity1 = velocities.col(index1);
            const Eigen::Vector3d& velocity2 = velocities.col(index2);
            const double damping1 = direction.dot(velocity1);
            const double damping2 = direction.dot(velocity2);

            const Eigen::Vector3d dampingForce1 = -dampingConstant * damping1 * direction;
            const Eigen::Vector3d dampingForce2 = -dampingConstant * damping2 * direction;

            forces.col(index1) += force + dampingForce1;
            forces.col(index2) += -force + dampingForce2 ;
        }

        // Adding Bending Springs
        for (geometrycentral::surface::Edge e : mesh.edges())
        {
            if (e.isBoundary())
            {
                continue;
            }

            geometrycentral::surface::Halfedge he = e.halfedge();
            geometrycentral::surface::Vertex vertex1 = he.next().next().vertex();
            geometrycentral::surface::Vertex vertex2 = he.twin().next().next().vertex();

            const size_t index1 = geometry.vertexIndices[vertex1];
            const size_t index2 = geometry.vertexIndices[vertex2];

            const Eigen::Vector3d diffVector = GCVec3toEigenVec3(geometry.vertexPositions[index2]) - GCVec3toEigenVec3(geometry.vertexPositions[index1]);
            const Eigen::Vector3d direction = diffVector.normalized();
            const double forceMagnitude = (diffVector.norm() - bendingLengths[e]) * bendingSpringConst;
            const Eigen::Vector3d force = forceMagnitude * direction;

            const Eigen::Vector3d& velocity1 = velocities.col(index1);
            const Eigen::Vector3d& velocity2 = velocities.col(index2);
            const double damping1 = direction.dot(velocity1);
            const double damping2 = direction.dot(velocity2);

            const Eigen::Vector3d dampingForce1 = -dampingConstant * damping1 * direction;
            const Eigen::Vector3d dampingForce2 = -dampingConstant * damping2 * direction;

            forces.col(index1) += force + dampingForce1;
            forces.col(index2) += -force + dampingForce2;
        }

        // Adding Gravity
        for (geometrycentral::surface::Vertex v : mesh.vertices())
        {
            const size_t index = geometry.vertexIndices[v];
            forces.col(index) +=  mass(index) * gravity;
        }

        // Explicit Euler
        auto massInverse = FlattenMatrixN3(invMass).asDiagonal();

        Eigen::Map<Eigen::VectorXd> f = FlattenMatrixN3(forces);

        FlattenMatrixN3(velocities) += timeStep * (massInverse * f);
        positions += timeStep * FlattenMatrixN3(velocities);

        // Saving to OBJ
        const std::string fileName = std::to_string(frame);
        if (!std::__fs::filesystem::exists("./output"))
        {
            std::__fs::filesystem::create_directory("./output");
        }
        geometrycentral::surface::writeSurfaceMesh(mesh, geometry,
                                                   "./output/" + fileName + ".obj");
    }

    std::cout<<"\n";
}


Eigen::Map<Eigen::Vector3d> GCVec3toEigenVec3(geometrycentral::Vector3& vector)
{
    return Eigen::Map<Eigen::Vector3d>(&vector.x, 3);
}

Eigen::Map<Eigen::VectorXd> FlattenMatrixN3(Eigen::Matrix<double, 3, Eigen::Dynamic>& vector)
{
    return Eigen::Map<Eigen::VectorXd>(vector.data(), 3 * vector.cols());
}

size_t pointToIndex(const geometrycentral::Vector3& point,
                  geometrycentral::surface::VertexPositionGeometry& geometry,
                  geometrycentral::surface::SurfaceMesh& mesh)
{
    size_t closestIndex = -1;
    double closestDistanceSq = std::numeric_limits<double>::max();
    for (geometrycentral::surface::Vertex v : mesh.vertices())
    {
        double distanceSq = (point - geometry.vertexPositions[v]).norm();
        if (distanceSq < closestDistanceSq)
        {
            closestIndex = geometry.vertexIndices[v];
            closestDistanceSq = distanceSq;
        }
    }
    return closestIndex;
}

void printProgress(double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}
