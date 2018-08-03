%% Function library
% Highly recommend using structures, e.g. a matlab structure "element"
% which describes the complete system where each object in the array of
% "element" is an element of the system, and implementing separate fields
% for each property of each object instance. Typically I was using six
% fields: element(n).k (the local stiffness matrix), element(n).L,
% element(n).theta, element(n).T, and element(n).u1 and element(n).u2 to
% record the ordinal number of the nodes at the ends of the nth element of
% the structure, e.g. if element 3 connected nodes 2 and 4, it would be
% element(3).u1 = 2, element(3).u2 = 4. 

%% Ordinary lirary functions
% zeros(), atan2(), norm(), and length().

%% Local matrices
% Used to build the local matrix for a bar element. This function can be
% adapted for springs by making area A = L = 1 and implementing the elastic
% modulus as hook's modulus. This function allows for arbitrary rotation.
function barElementMat = buildLocalBarMatrixAnyOrientation(L,angle,A,E)
    barElementMat = A*E/L*[localBarRotation(angle),-1*localBarRotation(angle); ... 
    -1*localBarRotation(angle),localBarRotation(angle)];
end

% Used to build the local matrix for a beam element. Rotation needs 
function beamElementMat = beamElement(E,I,L)
    beamMat = E*I/L^3*[12 6*L -12 6*L; 6*L 4*L^2 -6*L 2*L^2; -12 -6*L 12 -6*L; 6*L 2*L^2 -6*L 4*L^2];
end

% Used to build the local matrix for a frame element
function frameStiffnessMat = frameElement(A,E,L,I)
    c1 = A*E/L;
    c2 = E*I/L^3;
    frameStiffnessMat = [c1, 0, 0, -c1, 0, 0; 0, 12*c2, 6*c2*L, 0, -12*c2, 6*c2*L;...
        0, 6*c2*L, 4*c2*L^2, 0, -6*c2*L, 2*c2*L^2;...
        -c1, 0, 0, c1, 0, 0; 0, -12*c2, -6*c2*L, 0, 12*c2, -6*c2*L;...
        0, 6*c2*L, 2*c2*L^2, 0, -6*c2*L, 4*c2*L^2];
end

%% Superpositions

% Used to superposition either a bar or beam element. Here, matrix = local
% stiffness matrix, u1 is the lower ordinal number of the node, u2 is the
% higher ordinal number, dim is the degrees of freedom/rank of the matrix
% space. The higher dimensional space of frames (3 vs 2 for each node)
% required a second function (could've parameterized an argument but it was
% easier to copy-paste).
function localSuperMat = globalizeMatrix(matrix,u1,u2,dim)
    localSuperMat = zeros(dim,dim);
    localSuperMat(2*u1-1:2*u1,2*u1-1:2*u1) = matrix(1:2,1:2);
    localSuperMat(2*u2-1:2*u2,2*u1-1:2*u1) = matrix(3:4,1:2);
    localSuperMat(2*u1-1:2*u1,2*u2-1:2*u2) = matrix(1:2,3:4);
    localSuperMat(2*u2-1:2*u2,2*u2-1:2*u2) = matrix(3:4,3:4);
end

function localSuperMat = globalizeFrameMatrix(matrix,u1,u2,dim)
    localSuperMat = zeros(dim,dim);
    localSuperMat(3*u1-2:3*u1,3*u1-2:3*u1) = matrix(1:3,1:3);
    localSuperMat(3*u2-2:3*u2,3*u1-2:3*u1) = matrix(4:6,1:3);
    localSuperMat(3*u1-2:3*u1,3*u2-2:3*u2) = matrix(1:3,4:6);
    localSuperMat(3*u2-2:3*u2,3*u2-2:3*u2) = matrix(4:6,4:6);
end

%% Rotations

function rMat = frameRotationMat(theta)
    rMat = zeros(6,6);
    rMat(1:2,1:2) = rotationMat(theta);
    rMat(4:5,4:5) = rotationMat(theta);
    rMat(3,3) = 1;
    rMat(6,6) = 1;
%     matlab's evaluation of sin(0) != 0, I had to assert an
%     epsilonic/machine precision tolerance.
    machinePrecision = find(abs(rMat) < 1e-15);
    rMat(machinePrecision) = 0;
end

function localBarRotationMat = localBarRotation(theta)
    localBarRotationMat = rotationMat(theta)*rotationMat(theta)';
 
% matlab has problems with numerical precision for some of the products
% above at particular angles, need to assert tolerance/zero condition on
% points approaching machine precision.
    B = abs(localBarRotationMat) < 1e-15;
    localBarRotationMat(B) = 0;
end

function rMat = rotationMat(theta)
    rMat = [cos(theta),sin(theta);-sin(theta),cos(theta)];
end

%% Equivalent nodal loading
% These are the various formulas for calculating the equivalent nodal
% forces.

function array = distributedLoadCase1(L,P)
    array = zeros(4,1);
    array(1) = -P/2;
    array(2) = -P*L/8;
    array(3) = -P/2;
    array(4) = P*L/8;
end

function array = distributedLoadCase2(L,P,a,b)
    array = zeros(4,1);
    array(1) = -P*b^2*(L+2*a)/L^3;
    array(2) = -P*a*b^2/L^2;
    array(3) = -P*a^2*(L+2*b)/L^3;
    array(4) = P*a^2*b/L^2;
end

function array = distributedLoadCase3(L,P,alpha)
    array = zeros(4,1);
    array(1) = -P;
    array(2) = -P*L*alpha(1-alpha);;
    array(3) = -P;
    array(4) =P*L*alpha(1-alpha);
end

function array = distributedLoadCase4(w,L)
    array = zeros(4,1);
    array(1) = -w*L/2;
    array(2) = -w*L^2/12;
    array(3) = -w*L/2;
    array(4) = w*L^2/12;
end

function array = distributedLoadCase5(w,L)
    array = zeros(4,1);
    array(1) = -7*w*L/20;
    array(2) = -w*L^2/20;
    array(3) = -3*w*L/20;
    array(4) = w*L^2/30;
end

function array = distributedLoadCase6(w,L)
    array = zeros(4,1);
    array(1) = -w*L/4;
    array(2) = -5*w*L^2/96;
    array(3) = -w*L/4;
    array(4) = 5*w*L^2/96;
end

function array = distributedLoadCase7(w,L)
    array = zeros(4,1);
    array(1) = -13*w*L/32;
    array(2) = -11*w*L^2/192;
    array(3) = -3*w*L/32;
    array(4) = 5*w*L^2/192;
end

function array = distributedLoadCase8(w,L)
    array = zeros(4,1);
    array(1) = -w*L/3;
    array(2) = -w*L^2/15;
    array(3) = -w*L/3;
    array(4) = w*L^2/15;
end

function array = distributedLoadCase9(M,L,a,b)
    array = zeros(4,1);
    array(1) = -M(a^2+b^2-4*a*b-L^2)/L^3;
    array(2) = -M*b*(2*a-b)/L^2;
    array(3) = -M*(a^2+b^2-4*a*b-L^2)/L^3;
    array(4) = M*a*(2*b-a)/L^2;
end

