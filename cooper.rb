#!/usr/bin/env ruby

PAR_INPUT = "#{ENV['HOME']}/3dee_program/dd_version/forces/dp.dat"

PMASS = 938.2723

def main
  if ARGV.size != 3
    print_usage
    exit 1
  end

  a   = ARGV[0].to_i
  tp  = ARGV[1].to_f
  pot = ARGV[2]

  calc(a,tp,pot)
end

def print_usage
  puts "usage: ruby %s A E_p POT" % $0
  puts "POT: [EDAD1, EDAD2, EDAD3, EDAIC, EDAIO, EDAICA, EDAIZR, EDAIPB]"
end

def calc(a,e,pot)
  puts "parameters for A = %d, E = %.0f MeV, POT = %s" % [a,e,pot]

  af = a.to_f
  x  = 1000.0/calc_wp(a,e)
  y  = af/(af+20.0)

  pvec1 = [1.0, x, x**2, x**3, x**4, y, y**2, y**3]
  pvec2 = [x*y, x**2 * y, x * y**2]

  param = load_prm(pot)

  rfac = 1.0
  afac = 0.7

  rv1 = rfac * (inner_prod(param[1][0..7],pvec1) + inner_prod(param[13][0..2],pvec2))
  av1 = afac * (inner_prod(param[2][0..7],pvec1) + inner_prod(param[13][3..5],pvec2))
  puts "(R_VR, a_VR) = (%f, %f)"%[rv1,av1]

  rv2 = rfac * (inner_prod(param[4][0..7],pvec1) + inner_prod(param[14][0..2],pvec2))
  av2 = afac * (inner_prod(param[5][0..7],pvec1) + inner_prod(param[14][3..5],pvec2))
  puts "(R_VI, a_VI) = (%f, %f)"%[rv2,av2]

  #vc = (inner_prod(param[3][0..7],pvec1) + inner_prod(param[15][0..2],pvec2))
  #wc = (inner_prod(param[6][0..7],pvec1) + inner_prod(param[15][3..5],pvec2))
  #puts "(V   , W   ) = (%f, %f)"%[vc,wc]

  rs1 = rfac * (inner_prod(param[7][0..7],pvec1) + inner_prod(param[16][0..2],pvec2))
  as1 = afac * (inner_prod(param[8][0..7],pvec1) + inner_prod(param[16][3..5],pvec2))
  puts "(R_SR, a_SR) = (%f, %f)"%[rs1,as1]

  rs2 = rfac * (inner_prod(param[10][0..7],pvec1) + inner_prod(param[17][0..2],pvec2))
  as2 = afac * (inner_prod(param[11][0..7],pvec1) + inner_prod(param[17][3..5],pvec2))
  puts "(R_SI, a_SI) = (%f, %f)"%[rs2,as2]

  #vs = (inner_prod(param[9][0..7],pvec1) + inner_prod(param[18][0..2],pvec2))
  #ws = (inner_prod(param[12][0..7],pvec1) + inner_prod(param[18][3..5],pvec2))
  #puts "(V_s , W_s ) = (%f, %f)"%[vs,ws]

end

def load_prm (pot)
  require 'scanf'

  flag = false
  param_raw = []

  File.foreach(PAR_INPUT) do |line|
    flag = true if line[1..-3] == pot
    next unless flag
    param_raw += line.scanf("%f, %f, %f, %f,")
  end

  param = []

  # this code only work for 1.8.7 or later
  #while (param_raw.size > 0)
  #  param << param_raw.shift(8)
  #end

  length = param_raw.size/8
  for i in 0..(length-1)
    p = 8*(i)
    q = 8*(i+1)
    #puts "i = %d, (%d, %d)"%[i, p, q]
    param << param_raw[p..q]
    #param[i].each {|x| print x, ", "}
    #puts ""
    #puts "..........."
  end

  if param.length == 0
    puts "Error: Parameter Set named \'%s\' not found!" % pot
    exit 1
  end

  param
end

def calc_wp (a, tp)
  ep = tp + PMASS
  pp = Math::sqrt(tp*(tp+2*PMASS))
  beta = pp / (ep + a*PMASS)
  gamma = 1.0 / Math::sqrt((1+beta)*(1-beta))
  gamma * (ep - beta * pp)
end

def inner_prod (v1,v2)
  sum = 0.0
  #puts ".............."
  #puts v1
  #puts "-----"
  #puts v2
  #puts "========="
  v1.zip(v2).each do |e1,e2|
    sum += e1 * e2
  end
  sum
end

main
#end
