function dis =  calDistance(posi1,posi2)

point1 = locationToCartesian(posi1(1),posi1(2),0);
point2 = locationToCartesian(posi2(1),posi2(2),0);
dis = norm(point1-point2);
end 