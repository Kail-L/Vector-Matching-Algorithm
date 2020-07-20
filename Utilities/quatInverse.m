function [q_inv] = quatInverse(q)

% inverts a 1x4 quaternion (row vector)
%returns a 4x1 inverted quaternion (column vector)

q_a = [q(1,1),-q(1,2),-q(1,3),-q(1,4)]/norm(q);


q_inv = q_a';
