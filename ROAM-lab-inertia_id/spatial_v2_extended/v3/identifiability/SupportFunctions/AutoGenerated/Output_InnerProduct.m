function C1 = Output_InnerProduct(in1,in2)
%OUTPUT_INNERPRODUCT
%    C1 = OUTPUT_INNERPRODUCT(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    19-Dec-2017 16:00:33

v11 = in1(1,:);
v12 = in1(2,:);
v13 = in1(3,:);
v14 = in1(4,:);
v15 = in1(5,:);
v16 = in1(6,:);
v21 = in2(1,:);
v22 = in2(2,:);
v23 = in2(3,:);
v24 = in2(4,:);
v25 = in2(5,:);
v26 = in2(6,:);
C1 = [v14.*v24+v15.*v25+v16.*v26,-v12.*v26+v13.*v25+v15.*v23-v16.*v22,v11.*v26-v13.*v24-v14.*v23+v16.*v21,-v11.*v25+v12.*v24+v14.*v22-v15.*v21,v11.*v21,v12.*v22,v13.*v23,v12.*v23+v13.*v22,v11.*v23+v13.*v21,v11.*v22+v12.*v21];
