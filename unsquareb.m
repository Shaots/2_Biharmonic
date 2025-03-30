function [q,g,h,r]=unsquareb(p,e,u,time)
%UNSQUAREB	Boundary condition data.
%
%
%
bl=[
2 2 2 2
0 0 0 0
1 1 1 1
3 3 3 3
1 1 1 1
1 1 1 1
1 1 1 1
1 1 1 1
48 48 48 48
49 49 49 49
101 101 101 101
52 52 52 52
48 48 48 48
48 48 48 48
48 48 48 48
48 48 48 48
48 48 48 48
48 48 48 48
48 48 48 48
48 48 48 48
49 49 49 49
48 48 48 48
48 48 48 48
49 49 49 49
48 48 48 48
48 48 48 48
];

if any(size(u))
  [q,g,h,r]=pdeexpd(p,e,u,time,bl);
else
  [q,g,h,r]=pdeexpd(p,e,time,bl);
end
