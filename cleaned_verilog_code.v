// VerilogA for PLL_models, First_accumulator_new, veriloga
`include "constants.vams"
`include "disciplines.vams"
module First_accumulator_new(in,out,carry,clk);
inout [23:0] out;
inout carry,clk,in;
electrical [23:0] out;
electrical carry,clk,in;
parameter real VDD = 1;
integer old_value;
integer out_dec;
integer frac1;
integer sum ;
integer x;
integer out_temp[0:23];



real V_th;
integer quotient[0:23] ;
integer count;
integer count2;
integer count1;
integer temp5;
analog begin frac1 = V(in);
@(initial_step) old_value = 0;
V_th = VDD/2;
@(cross (V(clk) - V_th, 1)) begin
generate count2 (0,23,+1) begin
@(initial_step) V(out[count2]) <+ 0;
end
@(initial_step) V(out[23]) <+ 0; V(carry) <+ 0;
out_dec = frac1 + old_value;
old_value = 0;
generate count (0, 23, +1) begin
quotient [count] = out_dec/2;
if(out_dec % 2 != 0) out_temp[count] = 1;
else out_temp[count] = 0;
out_dec = quotient[count];
old_value = old_value + pow(2,count) * out_temp[count];


if (quotient[23] != 0 && quotient[23] % 2 == 1) temp5 = 1;
else temp5 = 0;
end
end
generate count1 (0,23,+1) begin V(out[count1]) <+ out_temp[count1]; end
V(carry) <+ temp5; end
endmodule


// VerilogA for PLL_models, Second_accumulator, veriloga
`include "constants.vams"
`include "disciplines.vams"
module Second_accumulator(in,out,carry,clk);
inout [23:0] out;
inout [23:0] in;
inout carry,clk;
electrical [23:0] out;
electrical [23:0] in;
electrical carry,clk;
parameter real VDD = 1.2;
integer old_value;


integer out_dec;
integer in_dec;
real V_th;
integer quotient[0:23] ;
integer count;
integer count1;
integer count2;
integer out_temp[0:23];
integer out_temp2;
analog begin
@(initial_step) old_value = 0;
@(initial_step) out_dec = 0;
@(initial_step) in_dec =0;
V_th = VDD/2;
@(cross (V(clk) - V_th, 1)) begin V(carry) <+ 0;
in_dec = 0;
generate count1 (0,23,+1) begin
in_dec = in_dec + pow(2,count1) * V(in[count1]); end
out_dec = in_dec + old_value; old_value = 0;
generate count (0,23,+1) begin quotient[count] = out_dec/2;


if(out_dec % 2 != 0) out_temp[count] = 1;
else out_temp[count] = 0;
out_dec = quotient[count];
old_value = old_value + pow(2,count) * out_temp[count];
if (quotient[23] != 0 && quotient[23] % 2 == 1) out_temp2 = 1;
else out_temp2 = 0; end
end
generate count2 (0,23,+1) begin V(out[count2]) <+ out_temp[count2];
end
V(carry) <+ out_temp2; end
Endmodule

// VerilogA for PLL_models, Third_accumulator, veriloga
`include "constants.vams"
`include "disciplines.vams"
module Third_accumulator(in,carry,clk);
inout [23:0] in;
inout carry,clk;
electrical [23:0] in;
electrical carry,clk;
parameter real VDD = 3.3;
integer old_value;
integer out_dec; integer in_dec; integer out[0:23]; real V_th;
integer quotient[0:23] ; integer count;
integer count1; integer out_temp2; analog begin
@(initial_step) old_value = 0; V_th = VDD/2;
@(cross (V(clk) - V_th, 1)) begin V(carry) <+ 0;
in_dec = 0;
generate count1 (0,23,+1) begin
in_dec = in_dec + pow(2,count1) * V(in[count1]); end
out_dec = in_dec + old_value; old_value = 0;
generate count (0,23,+1) begin quotient[count] = out_dec/2;
if(out_dec % 2 != 0) out[count] = 1;
else out[count] = 0; out_dec = quotient[count];
old_value = old_value + pow(2,count) * out[count];
if (quotient[23] != 0 && quotient[23] % 2 == 1) out_temp2 = 1;
else out_temp2 = 0; end
end
V(carry) <+ out_temp2; end
endmodule
// VerilogA for PLL_models, Delayy, veriloga

`include "constants.vams"
`include "disciplines.vams"
module Delayy(in,out,clk);
inout in,out,clk;
electrical in,out,clk;
parameter real VDD = 1.2;
real V_th;
integer old_value;
integer out_temp;
analog begin V_th = VDD / 2;
@(initial_step) old_value = 0; @(cross(V(clk) - V_th, 1)) begin
out_temp = old_value; old_value = V(in);
end
V(out) <+ out_temp; end
Endmodule

// VerilogA for PLL_models, first_adder_subtractor, veriloga

`include "constants.vams"
`include "disciplines.vams"
module first_adder_subtractor(in1,in2,out, clk);
inout in1,in2,out, clk;
electrical in1,in2,out, clk;
parameter real VDD = 1.2; real V_th;
integer temp; analog begin V_th = VDD /2;
@(cross(V(clk) - V_th, 1))
temp = V(in1) -V(in2); V(out) <+ temp;
end
endmodule

// VerilogA for PLL_models, first_adder, veriloga
`include "constants.vams"
`include "disciplines.vams"
module first_adder(in1,in2,out,clk);
inout in1,in2,out,clk;
electrical in1,in2,out,clk;
parameter real VDD = 1.2;
real V_th;
integer temp; analog begin V_th = VDD / 2;
@(cross(V(clk) - V_th, 1)) temp = V(in1) + V(in2);
V(out) <+ temp;
end endmodule
// VerilogA for PLL_models, first_adder_subtractor, veriloga

`include "constants.vams"
`include "disciplines.vams"
module first_adder_subtractor(in1,in2,out, clk);
inout in1,in2,out, clk;
electrical in1,in2,out, clk;
parameter real VDD = 1.2; real V_th;
integer temp; analog begin V_th = VDD /2;
@(cross(V(clk) - V_th, 1))
temp = V(in1) - V(in2); V(out) <+ temp;
end
endmodule

// VerilogA for PLL_models, final_adder, veriloga

`include "constants.vams"
`include "disciplines.vams"
module final_adder(in1,in2,clk,out);
inout in1, in2, clk;
inout [2:0] out;
electrical in1,in2, clk;
electrical [2:0] out;
parameter real VDD = 1.2;
integer sum;
analog begin

@(cross(V(clk) -(VDD/2),1)) begin sum = V(in1) + V(in2);
end
if (sum == -3) begin
V(out[0]) <+ 1;
V(out[1]) <+ 0;
V(out[2]) <+ 1;
end
else if (sum == -2) begin
V(out[0]) <+ 0;
V(out[1]) <+ 1;
V(out[2]) <+ 1;
end
else if (sum == -1) begin
V(out[0]) <+ 1;
V(out[1]) <+ 1;
V(out[2]) <+ 1;
end
else if (sum == 0) begin
V(out[0]) <+ 0;

V(out[1]) <+ 0;
V(out[2]) <+ 0;

end
else if (sum == 1) begin
V(out[0]) <+ 1;
V(out[1]) <+ 0;
V(out[2]) <+ 0;
end
else if (sum == 2) begin
V(out[0]) <+ 0;
V(out[1]) <+ 1;
V(out[2]) <+ 0;
end
else if (sum == 3) begin
V(out[0]) <+ 1;
V(out[1]) <+ 1;
V(out[2]) <+ 0;
end
else if (sum == 4) begin
V(out[0]) <+ 0;
V(out[1])<+0;
V(out[2])+1;
End
End
Endmodule


Code Verilog A của các khối FS, CCG, RAG và bộ tích luỹ



// VerilogA for PLL_models, Fraction_separation, veriloga
`include "constants.vams"
`include "disciplines.vams"

module Fraction_separation(Frac, N, f);
inout Frac, N, f;
electrical Frac, N, f;

real c;
real d;
integer a;

analog begin
c = V(Frac);
a = floor(c);
d = c - a;
V(N) <+ a;
V(f) <+ d;
end
endmodule


// VerilogA for PLL_models, control-code-generator, veriloga
`include "constants.vams"
`include "disciplines.vams"
module control_code_gen(f,i, k) ;
inout f,i,k;
electrical f, i, k; real b;
real e;
integer a;
real d;
integer m;
analog begin b = V(f);
e = b * 8;
a = floor(e); d = e - a; V(i) <+ a;
m = d * pow(2,24);
V(k) <+ floor(m); end
endmodule

// VerilogA for PLL_models, required-av-gen, veriloga

`include "constants.vams"
`include "disciplines.vams"
module required-av-gen(i, in ,out, B) ;
inout i, out, B;
inout [2:0] in;
electrical i, out, B;
electrical [2:0] in;
integer in_dec;
integer x;
integer y;
integer z; analog begin
if (V(in[0]) == 0 && V(in[1]) == 0 && V(in[2]) == 0)
in_dec = 0;
else if (V(in[0]) == 1 && V(in[1]) == 0 && V(in[2]) == 0) in_dec = 1;
else if (V(in[0]) == 0 && V(in[1]) == 1 && V(in[2]) == 0 ) in_dec = 2;
else if (V(in[0]) == 1 && V(in[1]) == 1 && V(in[2]) == 0 ) in_dec = 3;
else if (V(in[0]) == 0 && V(in[1]) == 0 && V(in[2]) == 1) in_dec = 4;

else if (V(in[0]) == 1 && V(in[1]) == 1 && V(in[2]) == 1) in_dec = -1;
else if (V(in[0]) == 0 && V(in[1]) == 1 && V(in[2]) == 1)

in_dec = -2;
else if (V(in[0]) == 1 && V(in[1]) == 0 && V(in[2]) == 1) in_dec = -3;
x= V(i);
y = x + in_dec; if(y <= 8) z = 0; if (y >8) begin z = 1;
y = y -8; end
V(out) <+ y; V(B) <+ z;
end
endmodule
// VerilogA for PLL_models, accumulator_2op, veriloga

`include "constants.vams"
`include "disciplines.vams"
module accumulator_2op(in,out1, out2,clk,C);
inout in;
inout clk;
inout [2:0] out1;
inout [2:0] out2;
electrical in;
electrical clk;
electrical [2:0] out1;
electrical [2:0] out2;
parameter real VDD = 1.5;
integer old_value;
integer old_value2;
integer in_dec;
integer out_dec;
integer count;
integer out_temp [0:2];
integer out_temp1[0:2];
integer quotient [0:2];
integer quotient1[0:2]; integer n;
integer r;
integer f;
integer M;
integer y;
integer S;
integer A;
integer count1;

analog begin
@(initial_step) old_value = 0; @(initial_step) r = 0; @(initial_step) in_dec = 0; @(initial_step) out_temp[0] = 0; @(initial_step) out_temp [1] = 0; @(initial_step) out_temp[2] = 0; @(initial_step) out_temp1[0] = 0; @(initial_step) out_temp1[1] = 0; @(initial_step) out_temp1[2] = 0; @(initial_step) S = 0; @(initial_step) Z = 0;
@(cross(V(clk) - (VDD/2) ,-1))
begin in_dec = V(in);
if(in_dec< 0) beginA = -1;in_dec = in_dec+8; end
if(in_dec <= 4)
out_dec = in_dec + old_value; if(in_dec > 4)
out_dec = 4 + old_value; if (out_dec == -3)
out_dec = 5;
if (out_dec == -2) out_dec= 6;
if (out_dec == -1) out_dec = 7;
if (out_dec >7) out_dec = out_dec - 8;
M = out_dec;
generate count (0,2,-1) begin
quotient[count] = out_dec/2;
if(out_dec % 2 != 0) out_temp[count] = 1;
else out_temp[count] = 0;
out_dec = quotient[count];
end
end
old_value = M;
V(out1[0]) <+ out_temp[0];
V(out1[1]) <+ out_temp[1];
V(out1[2]) <+ out_temp[2];
@(cross(V(clk) - (VDD/2) ,+1)) begin in_dec = V(in);
if (in_dec <= 4) begin
out_temp1 [0] = out_temp [0];
out_temp1 [1] = out_temp [1];
out_temp1 [2] = out_temp [2];
S = old_value; end
if(in_dec > 4) begin y = in_dec - 4;
r = y + old_value ; if (r == -3)
r = 5;
if (r == -2) r = 6;
if (r == -1) r = 7;
if (r >7) r = r - 8; S = r;
generate count1 (0,2,+1) begin quotient1[count1] = r/2;
if(r % 2 != 0) out_temp1[count1] = 1;
else out_temp1[count1] = 0;
r = quotient1[count1];
end
end
end
old_value = S;
V(out2[0]) <+ out_temp1[0];
V(out2[1]) <+ out_temp1[1];
V(out2[2])< + out_temp1[2];
End
Endmudule




Code Verilog A của bộ tổng hợp tần số số nguyên – N

// VerilogA for PLL_models, PFD_model, veriloga
`include "constants.vams"
`include "disciplines.vams"

module PFD_model(ref_clk, fb_clk, up_out, down_out);

// *****************************
// **** Definitions *********

// define which signal is input, output and inout
inout ref_clk;
inout fb_clk;
inout up_out;
inout down_out;

// define that the signals are electrical
electrical ref_clk, fb_clk, up_out, down_out;

// define parameters to enter
// syntax:
// parameter type name = default value;
parameter real VDD = 1;
parameter real ttol = 10p; real td ;
parameter real tr = 140f; parameter real tf = 140f;
// define internal variables
integer state; // state=1 for down, -1 for up
// ***************************
//***** functionality*********
analog begin
// cross (condition, direction, tol) direction +1:rise, -1:fall, 0:
both
@(initial_step) td = 2n;
@(cross( V(ref_clk)-VDD/2,1,ttol))state=state-1;
@(cross( V(fb_clk)-VDD/2,1,ttol)) state=state+1;
if (state >1) state=1;
if (state <-1) state=-1;


V(down_out) <+ transition((state+1)/2*VDD,td,tr,tf); V(up_out) <+ transition((state-1)/2*VDD+VDD,td,tr,tf);
end // analog
endmodule
// VerilogA for PLL_models, CP_model, veriloga

`include "constants.vams"
`include "disciplines.vams"
module CP_model(Iout , Down , up);
inout Iout;
inout Down;
inout up;
electrical Iout, Down, up;
parameter real Ip = 20e-6;
parameter real V_max = 1;
parameter real V_min = 0;
parameter real Mis = 1e-9;
parameter real VDD = 1;
parameter real Delay = 1f;
parameter real tr= 1f;
parameter real tf = 1f; integer state;
real v_th; analog begin v_th = VDD/2;

@(cross( V(up) - v_th,-1)) state= -1;
@(cross( V(up) - v_th,1)) state = 0;
@(cross( V(Down) - v_th,-1)) state = 1;
@(cross(V(Down) - v_th, 1)) state = 0;
@(cross(V(Iout) - V_max,1)) state=0;
if (state>-1) state=0;
if (state <1) state=0;
I(Iout) <+ transition(Ip*(state)*(1+state*Mis),tr,tf);
end
endmodule
// VerilogA for PLL_models, vco_multiplephases_model, veriloga

`include "constants.vams"
`include "disciplines.vams"
module vco_multiplephases_model(in,out1,out2,out3,out4,out5,out6,out7,out8);
inout in,out1,out2,out3,out4,out5,out6,out7,out8;
electrical in,out1,out2,out3,out4,out5,out6,out7,out8;
//parameter real VDD = 3.3;
parameter real amp = 0.6;
parameter real offset = 0;
parameter real KVCO = 90e6;
parameter real fnom1 = 1.94e9;
parameter real start_time = 200u;
real freq;
real m; analog begin m = fnom1;
freq= m +KVCO*V(in);
V(out4)<+amp*sin(2 * `M_PI *freq *$abstime)+offset;

V(out5)<+amp*sin((2 * `M_PI *freq *$abstime)-(2*`M_PI/8))+offset;
V(out6)<+amp*sin((2 * `M_PI *freq *$abstime)-(2*2*`M_PI/8))+offset; V(out7)<+amp*sin((2 * `M_PI *freq *$abstime)-(3*2*`M_PI/8))+offset; V(out8)<+amp*sin((2 * `M_PI *freq *$abstime)-(4*2*`M_PI/8))+offset; V(out1)<+amp*sin((2 * `M_PI *freq *$abstime)-(5*2*`M_PI/8))+offset; V(out2)<+amp*sin((2 * `M_PI *freq *$abstime)-(6*2*`M_PI/8))+offset; V(out3)<+amp*sin((2 * `M_PI *freq *$abstime)-(7*2*`M_PI/8))+offset;
//The following part is made active if it is needed to change the tuning curve during simulation.
//@(timer(start_time))
//m = fnom2 ;
$bound_step(0.05 / freq); end
endmodule
// VerilogA for PLL_models, Divder_modified, veriloga
`include "constants.vams"
`include "disciplines.vams"
module Divder_modified(clock_out,clock_in, N, B, C);

inout clock_out, clock_in, N,B, C;
electrical clock_out,clock_in, N, B, C;
parameter real V_high= 1.5;
parameter real V_low =0;
parameter real td=15f;


parameter real tr= 15f;
parameter real tf = 15f;
integer count;
integer n; real m;
analog begin
@(cross((V(clock_in)-0.3),1)) begin m= V(N) + V(B)+ V(C);
count = count + 1;
if(count>= m) count = 0 ;
n = (2*count >= m);
end
V(clock_out)<+transition(n? V_high:V_low,td, tr, tf);
end
endmodule
// VerilogA for PLL_models, Mux_2sel, veriloga
`include "constants.vams"
`include "disciplines.vams"
module Mux_2sel(sel1,sel2, clk,out1,out2,out3,out4,out5,out6,out7,out8,Y);
Pinout out1,out2,out3,out4,out5,out6,out7,out8,Y;
inout [2:0] sel1;
inout [2:0] sel2; inout clk;

electrical out1,out2,out3,out4,out5,out6,out7,out8,Y;
electrical [2:0] sel1;
electrical [2:0] sel2;
electrical clk;
parameter real VDD = 1;
integer sel_value1;
integer sel_value2;
integer sel_value;
analog begin
@(cross(V(clk)-(VDD/2), 1)) begin
sel_value1 = 1*V(sel1[0])+2*V(sel1[1])+4*V(sel1[2]);
sel_value = sel_value1;
end
@(cross(V(clk)-(VDD/2), -1)) begin
sel_value2 = 1*V(sel2[0])+2*V(sel2[1])+4*V(sel2[2]);
sel_value = sel_value2;
end
if (sel_value == 0) V(Y)<+ V(out4);
if (sel_value == 1) V(Y) <+ V(out5);
if (sel_value == 2) V(Y) <+ V(out6);
if (sel_value == 3) V(Y) <+ V(out7);
if (sel_value == 4) V(Y) <+ V(out8);
if (sel_value == 5) V(Y) <+ V(out1);
if (sel_value == 6) V(Y) <+ V(out2);
if (sel_value == 7) V(Y) <+ V(out3);
end
endmodule


