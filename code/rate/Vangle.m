function a=Vangle(u,v) 
% a=atan2(norm(cross(u,v)),dot(u,v));
a=acosd(dot(u,v)/(norm(u)*norm(v)));