% Loads the topcat library into matlab

clc; clear all; close all;
clear import;

javaaddpath('../target/topcat-1.0-SNAPSHOT.jar');
import topcat.*;

addpath('./examples');