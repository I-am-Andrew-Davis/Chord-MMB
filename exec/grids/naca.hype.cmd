# Mimic NASA turbulent flow over an airfoil reference
#
# grid size
$nRadial=32;
$nTheta=64;
# expected AMR (overestimaes are good)
$ref=16;
# solve number of lines on finest level
$nfTheta = $ref*$nTheta + 1;
$nfRadial = $ref*$nRadial;
#
create mappings
  # make the NACA airfoil (curve)
  Airfoil
    airfoil type
      naca
    lines
     $nfTheta*10
    mappingName
      airfoil-curve
  #pause
    exit
  # A nurbs formulation -- faster? and uniform point distribution
  nurbs (curve)
    interpolate from a mapping
      airfoil-curve
    lines
      $nfTheta*10
    mappingName
      airfoil-nurbs
    exit
#
  hyperbolic
#
    BC: left trailing edge
    BC: right trailing edge
    BC: bottom match to a mapping
      airfoil-nurbs
#
    points on initial curve $nfTheta
    lines to march $nfRadial
    distance to march 500
    marching options...
      volume smooths 500
      uniform dissipation coefficient .5
    geometric stretch factor 1.02
    generate
#
    smoothing...
    GSM:BC: left periodic
    GSM:BC: right periodic
    GSM:BC: top slide
    GSM:BC: bottom slide
    GSM:number of iterations 400
    GSM:number of elliptic smooths 10
    GSM:number of control function smooths 3
    GSM:smooth grid
#    lines
#      $nfRadial $nfTheta
#
#    plot ghost lines
#      0
#  pause
    mappingName
     airfoil
    boundary conditions
      -1 -1 1 2
    exit
# convert to a nurbs
#  nurbs (curve)
#    interpolate from mapping with options
#      airfoil
#      parameterize by index (uniform)
#      choose degree
#        5
#      done
#    mappingName
#      airfoil-nurbs-grid
#    exit
#
# Grid stretching towards the back of the airfoil
  stretch coordinates
    transform which mapping?
      airfoil
    Stretch r1:itanh  
    STP:stretch r1 itanh: layer 1 1 15 0
    Stretch r1:itanh   
    STP:stretch r1 itanh: layer 2 1 10 0.5
    mappingName
      airfoil-stretch
# pause
    exit 
#
# pause
 exit
#
# make an overlapping grid
#
generate an overlapping grid
    airfoil-stretch
#    airfoil
  done
  change parameters
    ghost points
      all
      2 2 2 2 2 2
  exit
  # pause
  compute overlap
  # pause
exit
#
save an overlapping grid
naca.hype.hdf
naca
exit
