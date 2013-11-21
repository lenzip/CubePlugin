CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(cuboPhysicsList)

FIND_PACKAGE(GGS REQUIRED)

ADD_LIBRARY(cuboPhysicsList SHARED src/cuboPhysicsList.cc src/cuboPhysicsListMessenger.cc)
