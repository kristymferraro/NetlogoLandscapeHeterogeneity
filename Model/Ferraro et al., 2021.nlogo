;running this model 15 years
;N is delayed for caracss and fecal input
;N is binned in groups of 50
;here, I will attempt to create a version of the model where there is only consumption, defication, and collective
;20% of carcass N ends up in patch
;set up sex ratio
;only 90% of females are pregnant

globals [total-dead adult-dead calf-dead N-Adult-dead-count N-calf-dead-count Total-Caribou Days Hour Day TEST  agedeath ] ; these are varibales that will be spit out in data sheet

breed [caribous caribou]

patches-own [
 patch-n ; nitrogen in a given patch
 alive? ; does the patch have enough nitrogen to have food?
  patch-counter-excretion ;(unsure if I need this - but maybe for diffusion or assimilation lag time)
  deposited-excretion-n ;the pool to keep the N that is deposted
  patch-counter-carcass
  deposited-carcass-n ;the pool to keep the N that is deposted
  Total-N-Depositied ; this pool will add up all the N that is deposited in the patch
  Total-N-Consumed
]

turtles-own [
 sex
mom?
 body-n ; amount of N in the body
 age ; age of the cairbou
 caribou-tick
 caribou-day ; used for aging, at day 90, cairbou age
 has-calf? ; does the caribou have a calf?
 is-calf? ;is the individual a calf?
 orient
 step
 daily-n
 yesterday-n
  excretion-n
  yesterday-n-left
 leader? ; identifying the random leaders
 closest-leader ; identifying the closest leader
 caribou-year
]



to setup
clear-all
set adult-dead  0
set calf-dead  0
set Total-Caribou 0
  ifelse Heterogeneity[
    ask patches [
      set patch-n random 1000 ;gives each patch a random amount of nitrogen between 0-1000
      set pcolor scale-color green patch-n 0 1000 ; set pcolor shade of green based on n (makes model run slower)
       if patch-n >= 0 [
        set pcolor brown ;the color is set to brown
        set alive? False ];and the patch is set to be dead
       if patch-n > 0 [ ;if the patch has > 0 nitrogen
        set pcolor scale-color green patch-n 0 1000 ;the patch is assigned a color
        set alive? True ];and set to be alive
      set patch-counter-excretion 0
      set patch-counter-carcass 0
  ]]
  [ask patches [
      set patch-n 500
      set pcolor scale-color green patch-n 0 1000
      set patch-counter-excretion 0
      set patch-counter-carcass 0
  ]]
create-caribous Population
  ask caribous [
    setxy random-xcor random-ycor
    set color red
    set size 2
    set age random 7 + 1
       if age = 1 [ set body-n random-normal 3450 3 ]
       if age = 2 [ set body-n random-normal 3750 3 ]
       if age > 2 [ set body-n random-normal 4000 3 ]
    set has-calf? false
    set is-calf? false
    set caribou-day 0
    set daily-n 0
    set leader? false
    set mom? false
    set caribou-year 0 ]
 ask n-of (count caribous * 0.7) caribous [
      set sex "female" ]
 ask caribous with [ sex != "female" ] [
       set sex "male" ]
reset-ticks

end


to go
  choose-leaders
   choose-moms
  ask turtles [
  move
  consume-n
  deficate
  end-day
  caribou-time
  count-caribous
  set caribou-tick caribou-tick + 1    ; each tick is 1/5 a day
  ]
ask patches [
  set-patch-N
  defecate-delay
  carcass-delay
  ]

tick
    if ticks >= 131400 [stop] ;15 years with 24 ticks a day
end


; = = = = = = = = == = = = = = = = = = = = Submodels for Both  = = = = = = = = = = = = = = = = = = = = = =
; TIME: Each day has 20 ticks to it.
to caribou-time
  if caribou-tick = 24 [
  set caribou-day caribou-day + 1
  set caribou-tick 0 ]
end

to end-day
  if caribou-tick = 24 [
  population-control
  N-related-death
  reproduce-caribou
  end-year
  set body-n body-n + (daily-n * .38) ; add the N that has been eated all day (35% of N eaten is added) (which is 14g if they are on average consuming 40g a day)
  set body-n body-n - 14 ; expected average daily N needs so subtract this for homeostasis
    set excretion-n (daily-n * .62) + 14 ;urine n is set to what is not digested (60%) plus what is used for metabolims
    set yesterday-n 0
    set yesterday-n excretion-n ;move what should be urinated to be urinated tomorrow
  set daily-n 0
    set yesterday-n-left yesterday-n
  ]
end

to end-year
  if caribou-day = 365 [ ;at the end of each year
  set age age + 1 ;each caribou gains an age
  set caribou-day 0 ;the day is set to 0
  set has-calf? false ;they no longer are bound to nursery groups
  set is-calf? false  ;all calves become adults
  set mom? false
  set caribou-year caribou-year + 1
  ]
end

to move ; using the Levy Walk
  ifelse Social-groups[ ;if the social groups switch is on then social groups will form and other will use levy-walk
  ifelse caribou-day >= 90 and caribou-day <= days-of-social + 90 [
   if leader? = true [ ;if a leader, they will use levy walk, unless on a dead pach, in which then they face a alive patch and move 1
      ifelse alive? = True [ levy-walk]
        [face  one-of patches with [ patch-n > 0]
          fd 1]]
   if has-calf? or is-calf? = true and leader? = false [ ;if not a leader, and a calf or mom, they follow leader
      follow-leader]
   if has-calf? or is-calf? = false [ ;if not a calf or mom, they use levy walk
        levy-walk]]
    [levy-walk]]
  [levy-walk] ; if social groups if off, we use levy walk
end

to levy-walk
   set orient turning-angle-range 0 360
    set heading heading + one-of [-1 1] * orient
    set step power-law-dist 1 min-move-length scale-exp
    fd step
end

to-report turning-angle-range [ #min #max ]  ; function for determining the initial turning angle preceding any given move as a random floating number within a specified range
  report #min + random-float ( #max - #min ) ; calculates and reports range of turning angles for turtles to undertake during movement
end ; end of turning angle range calculation and reporting procedure

to-report power-law-dist [ #norm-const #min-step-length #scale-exp ] ; function defining a power law distribution for LW based on the normalization constant, the
  ; minimum step length and the scaling exponent of the distribution
  set scale-exp #scale-exp ; sets the scaling exponent of the power law distribution function to a randomly generated scaling exponent variable
  let randomizer random-float 1
  let min-step-length ( #norm-const * #min-step-length * ( 1 - randomizer ) ^ ( - 1 / #scale-exp )) ; calculates the minimum step length executed by LW subsidies from
  ;the normalization constant, scaling exponent and minimum step length of the power law function
  report min-step-length ; produces the step length for the LW
end ; end of power-law distribution for calculation and reporting procedure for defining turtle LW step length

to choose-leaders
  if ticks = 2185 or ticks = 10921 or ticks = 19681 or ticks = 28441 or ticks = 37201 or ticks = 45961
  or ticks = 54721 or ticks = 63481 or ticks = 72241 or ticks = 81001 or ticks = 89761 or ticks = 98521
  or ticks = 107281 or ticks = 116041 or ticks = 124801 or ticks = 133561[ ; needs to be the tick after caribou have had calves
   ask n-of ((count caribous with [has-calf? = true]) * 0.5) caribous with [sex = "female"]
    [set leader? TRUE
      set color 25]
  ]
end


to follow-leader ;;
  if (count caribous with [leader? = true]) = 0 [
        ask one-of caribous with [leader? = false]
        [set leader? TRUE
        set color 25]]
  if has-calf? = true or is-calf? = true [
    let nearby-leaders turtles with [leader? = true] in-radius 5 ;; find nearby leaders within a radius of 5 patches
    set closest-leader min-one-of other turtles with [leader? = true] [distance myself] ; sets the closest learder to be which ever leader is closest
    ifelse any? nearby-leaders [ ; if there are any leaders in 10 patches use levy walk
    levy-walk]
    [move-to closest-leader]; if not, then move to the closet leader
  ]
end

to remove-social-groups
  if caribou-day = days-of-social + 90 [
  set has-calf? false ;they no longer are bound to nursery groups
  set is-calf? false  ;all calves become adults
  set leader? FALSE ; clear all previous leaders
  ]
end


to choose-moms
  if ticks = 1 or ticks = 8900 or ticks = 17700 or ticks = 26400 or ticks = 35200
  or ticks = 44000 or ticks = 52800 or ticks = 61500 or ticks = 70300 or ticks = 79000
  or ticks = 87800 or ticks = 96600 or ticks = 105300 or ticks = 114000 or ticks = 122941 [ ; needs to be the tick after caribou have had calves
    ask n-of ((count caribous with [sex = "female" and  age >= 3]) * 0.9) caribous with [sex = "female" and  age >= 3]
    [set mom? TRUE]
 ]
end



to reproduce-caribou  ;caribou procedure
 if caribou-day = 90 and mom? = true[
    set body-n body-n - 240 ;give away N that is in calf
   hatch 1 [ rt random-float 360 ;hatch an offspring
      fd 1 ;move it forward one
      set age 0 ;give it an age of 0
      set body-n 240
      set is-calf? true
      set mom? false
       set has-calf? true ;
   let male-female random-float 1
     ifelse  male-female < 0.5 [ set sex "female"] [ set sex "male"]
    ]
      set has-calf? true ;give adult a the calf
    ]
end

to consume-n  ; caribou procedure
  if patch-n <= 0 [ ;if the patch has no N
    set pcolor brown ;the color is brown
    set alive? False ];the patch is dead
    if alive? = True  [ ;older caribou end up on a rich patch
      if is-calf? = false [
        if patch-n >= 1.67 [ ;1.67
          set patch-n patch-n - 1.67 ;patch losses n
          set daily-n daily-n + 1.67
          set Total-N-Consumed Total-N-Consumed + 1.67 ] ;caribou eats lots!
        if patch-n <= 1.67 and patch-n > 0 [
          let left-n  patch-n
          set patch-n patch-n - left-n ;patch losses n
          set daily-n daily-n + left-n
          set Total-N-Consumed Total-N-Consumed + left-n  ] ;caribou eats only a little
      ]
      if is-calf? = true[
        if patch-n >= 2.92 [
          set patch-n patch-n - 2.92 ;patch losses n
          set daily-n daily-n + 2.92
          set Total-N-Consumed Total-N-Consumed + 2.92] ;caribou eats lots!
        if patch-n <= 2.92 and patch-n > 0 [
          let left-n  patch-n
          set patch-n patch-n - left-n ;patch losses n
          set daily-n daily-n + left-n
          set Total-N-Consumed Total-N-Consumed + left-n  ] ;caribou eats only a little
    ]]
end

to set-patch-N ; patch procedure
 if patch-n <= 0 [
    set pcolor brown ;the color is set to brown
    set alive? False ;and the patch is set to be dead
  ]
  if patch-n > 0 [ ;if the patch has > 0 nitrogen
    set pcolor scale-color green patch-n 0 1000 ;the patch is assigned a color
    set alive? True ;and set to be alive
  ]
end

to N-related-death
  if  age >= 2 [
  if body-n <= 3700 [ ;asking the caribou n poor to die
    let probability-survive .20 ; the prob of for the survival of the whole group  (here the death rate is 80%)
    let math 1 - probability-survive ^ (1 / 365)  ; the prob that any indiviual dies on any given day part
    let probability-die 1 / math ;  the prob that any indiviual dies on any given day part
    let death-random random probability-die ; setting it so that each indiviual has that probability of dying
     if death-random = 1 [
    set N-Adult-dead-count N-Adult-dead-count + 1
    set total-dead total-dead + 1
      deposit-body-n  ;this will only happen if cumlative is on
     set TEST body-n
        set agedeath age

    die
  ]]]
   if age = 1 [
  if body-n <= 3200 [ ;asking the caribou n poor to die
    let probability-survive .20 ; the prob of for the survival of the whole group  (here the death rate is 80%)
    let math 1 - probability-survive ^ (1 / 365)  ; the prob that any indiviual dies on any given day part
    let probability-die 1 / math ;  the prob that any indiviual dies on any given day part
    let death-random random probability-die ; setting it so that each indiviual has that probability of dying
     if death-random = 1 [
    set N-Adult-dead-count N-Adult-dead-count + 1
    set total-dead total-dead + 1
        deposit-body-n ;this will only happen if cumlative is on
    die
  ]]]

end

to population-control
if is-calf? = true [ ;if a calf
    let probability-survive .35 ; the prob of for the survival of the whole group  (here the death rate is 65%)
    let math 1 - probability-survive ^ (1 / 365)  ; the prob that any indiviual dies on any given day part
    let probability-die 1 / math ;  the prob that any indiviual dies on any given day part
    let death-random random probability-die ; setting it so that each indiviual has that probability of dying
     if death-random = 1 [
    deposit-body-n ;this will only happen if cumlative is on
    set total-dead total-dead + 1
    set calf-dead calf-dead + 1
    die]]
if age >= 1 and age <= 5  and sex = "female" [
    let probability-survive .80 ; the prob of for the survival of the whole group  (here the death rate is 20%)80
    let math 1 - probability-survive ^ (1 / 365)  ; the prob that any indiviual dies on any given day part
    let probability-die 1 / math ;  the prob that any indiviual dies on any given day part
    let death-random random probability-die ; setting it so that each indiviual has that probability of dying
     if death-random = 1 [
    deposit-body-n ;this will only happen if cumlative is on
    set total-dead total-dead + 1
    set adult-dead adult-dead + 1
    die ]]
if age >= 1 and age <= 5 and sex = "male" [
    let probability-survive .65 ; the prob of for the survival of the whole group  (here the death rate is 20%)80
    let math 1 - probability-survive ^ (1 / 365)  ; the prob that any indiviual dies on any given day part
    let probability-die 1 / math ;  the prob that any indiviual dies on any given day part
    let death-random random probability-die ; setting it so that each indiviual has that probability of dying
     if death-random = 1 [
    deposit-body-n ;this will only happen if cumlative is on
    set total-dead total-dead + 1
    set adult-dead adult-dead + 1
    die ]]
if age > 5 and age >= 11 and sex = "female" [
    let probability-survive .90 ; the probability of for the survival of the whole group  (here the death rate is 20%)
    let math 1 - probability-survive ^ (1 / 365)
    let probability-die 1 / math
    let death-random random probability-die
     if death-random = 1 [
     deposit-body-n ;this will only happen if cumlative is on
    set total-dead total-dead + 1
    set adult-dead adult-dead + 1
    die]]

if age > 5 and age >= 11 and sex = "male" [
    let probability-survive .75 ; the probability of for the survival of the whole group  (here the death rate is 20%)
    let math 1 - probability-survive ^ (1 / 365)
    let probability-die 1 / math
    let death-random random probability-die
     if death-random = 1 [
     deposit-body-n ;this will only happen if cumlative is on
    set total-dead total-dead + 1
    set adult-dead adult-dead + 1
    die]]
if age > 11[
    let probability-survive .10 ; the probability of for the survival of the whole group  (here the death rate is 20%)
    let math 1 - probability-survive ^ (1 / 365)
    let probability-die 1 / math
    let death-random random probability-die
     if death-random = 1 [
  deposit-body-n ;this will only happen if cumlative is on
      set total-dead total-dead + 1
      set adult-dead adult-dead + 1
      die]  ]
if age > 12[
  deposit-body-n ;this will only happen if cumlative is on
      set total-dead total-dead + 1
      set adult-dead adult-dead + 1
    die]

;if Population = 8 [
;  if count turtles > 12 [
;    let x count turtles - 10
;    ask n-of x turtles [
;      deposit-body-n ;this will only happen if cumlative is on
;      set total-dead total-dead + 1
;        if age  >= 1 [ set adult-dead adult-dead + 1]
;        if age  < 1 [ set calf-dead calf-dead + 1]
;       die]]]
;if Population = 72 [
;  if count turtles > 100 [
;    let x count turtles - 100
;    ask n-of x turtles [
;      deposit-body-n ;this will only happen if cumlative is on
;      set total-dead total-dead + 1
;      if age  >= 1 [ set adult-dead adult-dead + 1]
;      if age  < 1 [ set calf-dead calf-dead + 1]
;      die]]]
;if Population = 180 [
;  if count turtles > 252 [
;    let x count turtles - 252
;    ask n-of x turtles [
;      deposit-body-n ;this will only happen if cumlative is on
;      set total-dead total-dead + 1
;      if age  >= 1 [ set adult-dead adult-dead + 1]
;        if age  < 1 [ set calf-dead calf-dead + 1]
;      die]]]
end

to count-caribous
    set Total-Caribou count turtles
end

; = = = = = = = = == = = = = = = = = = = = Cumulative  = = = = = = = = = = = = = = = = = = = = = =
to deficate
 ifelse Interactions = "Cumulative"[ ;only do this if cumlative is on
let excretion-rate yesterday-n / 10 ;the caribou deficate 10 times a day from the pool the day before
set yesterday-n-left yesterday-n-left
  if ticks mod 2 = 0 [ ; every other tick
    let n-left yesterday-n-left - excretion-rate
 set yesterday-n-left n-left
    ask patch-here [ ; if the patch has an agent on it
    set deposited-excretion-n deposited-excretion-n + excretion-rate
    set Total-N-Depositied Total-N-Depositied + excretion-rate
  ]]]
  []
  end

to defecate-delay
  if  deposited-excretion-n > 0 [
  if ticks mod 24 = 0 [ ; every 20 ticks
      set patch-counter-excretion patch-counter-excretion + 1]]
    if patch-counter-excretion = 30  [
      set patch-n patch-n + deposited-excretion-n
      set patch-counter-excretion 0
      set deposited-excretion-n 0 ]
end


to deposit-body-n ;deposition of carcass N
  ifelse Interactions = "Cumulative"[ ;only do this if cumlative is on
  if is-calf? = false [
    let carcass-n body-n
    ask patch-here [
     set deposited-carcass-n (carcass-n * .2)
      set Total-N-Depositied Total-N-Depositied + (carcass-n * .2)] ;adding 20% of nitrogen from body to the patch (total N in body)
    ]
  if is-calf? = True [
    let carcass-n body-n
    ask patch-here [
    set deposited-carcass-n  (carcass-n * .2)] ;adding nitrogen to the patch (total N in body)
      set Total-N-Depositied Total-N-Depositied + (carcass-n * .2) ]]
    []
end


to carcass-delay
  if  deposited-carcass-n  > 0 [
  if ticks mod 24 = 0 [ ; every 20 ticks
      set  patch-counter-carcass  patch-counter-carcass + 1]]
    if  patch-counter-carcass = 120  [
      set patch-n patch-n + deposited-carcass-n
      set patch-counter-carcass 0
      set deposited-carcass-n  0]
     ifelse patch-n > 0 [ ;bring back the patches to life if they get more N
        set pcolor scale-color green patch-n 1 1000
        set alive? TRUE ]
      [set pcolor brown
      set alive? FALSE ]
end
@#$#@#$#@
GRAPHICS-WINDOW
288
14
875
602
-1
-1
5.733
1
10
1
1
1
0
1
1
1
-50
50
-50
50
1
1
1
ticks
30.0

SLIDER
2
241
174
274
Population
Population
0
180
50.0
1
1
NIL
HORIZONTAL

BUTTON
6
10
72
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
8
51
71
84
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
84
10
163
55
dead-count
total-dead
17
1
11

MONITOR
80
60
163
105
NIL
Total-Caribou
17
1
11

SWITCH
8
206
152
239
Social-groups
Social-groups
1
1
-1000

SLIDER
4
289
176
322
scale-exp
scale-exp
0
3
2.0
1
1
NIL
HORIZONTAL

SLIDER
2
327
174
360
min-move-length
min-move-length
0
3
3.0
1
1
NIL
HORIZONTAL

SLIDER
2
365
174
398
days-of-social
days-of-social
0
180
180.0
90
1
NIL
HORIZONTAL

MONITOR
174
26
231
71
Days
[caribou-day] of one-of caribous
17
1
11

PLOT
1102
14
1302
145
Population
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"pen-1" 1.0 0 -7500403 true "" "plot count turtles"

PLOT
883
13
1083
163
Number of Dead Patches
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count patches with [alive? = false]"
"pen-1" 1.0 0 -7500403 true "" "plot count patches with [patch-n <= 0]"

PLOT
889
180
1089
330
Adult Average daily-N
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [daily-n] of caribous with [is-calf? = false]"

PLOT
898
351
1098
501
Calf daily N
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [daily-n] of caribous with [is-calf? = true]"

PLOT
1110
340
1310
490
Calf Body N Ave
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [body-n] of caribous with [is-calf? = TRUE]"

PLOT
1105
171
1305
321
Adult Average Body-N
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [body-n] of caribous with [is-calf? = FALSE]"

MONITOR
174
77
231
122
Year
[caribou-year] of one-of caribous
17
1
11

SWITCH
7
167
150
200
Heterogeneity
Heterogeneity
0
1
-1000

MONITOR
20
411
162
456
NIL
N-Adult-dead-count
17
1
11

MONITOR
21
462
152
507
NIL
N-calf-dead-count
17
1
11

CHOOSER
10
119
148
164
Interactions
Interactions
"Consumption" "Cumulative"
1

@#$#@#$#@
## WHAT IS IT?

This model seeks to explore how caribou move nutrients around in the summer range considering both the grazing and release of nutreints through waste and carcass depostion. 

## HOW IT WORKS

We constructed a spatially explicit individual-based model in NetLogo version 6.1.1. We designed an environment to represent a 600 × 600 square grid, wrapped vertically and horizontally. Each patch (10 m2 size) represented the nitrogen within vegetation that could be consumed and as well as the nitrogen deposited by agents that then is incorporated into the soil, and subsequently the plant matter. 

Caribou – The main agent of this model simulates the characteristics of Rangifer tarandus (caribou). In the this cumulative versions of the model, the caribou removes nitrogen from each patch it visits by foraging on the vegetation. If a caribou visits a patch with available “plant” nitrogen, the caribou gains nitrogen; however, if a caribou visits a patch without available “plant” nitrogen, the caribou gains no nitrogen. As well, if a caribou visits a patch without nitrogen, in their next step they will reorient themselves towards a patch with plant nitrogen. At the end of each day, caribou retain a fixed percentage of the nitrogen consumed throughout the day. Additionally, nitrogen needed for metabolic process is moved from the body pool and is excreted the following day. The cumulative model also seeks to also explore the return of nutrient to the system by animals. Therefore, in the cumulative model, a fixed percentage of the daily consumed nitrogen is also set aside to be deposited the following day, simulating excretion Additionally, at death, the caribou carcass will deposit the nitrogen within their body on the soil patch in which they die

Patch environment – The patches, i.e. the model environment, represent areas that contain nitrogen in plant tissue. In the consumption model caribou can consume the nitrogen within a given patch while in the cumulative model the caribou can both consume and excrete nitrogen within a given patch. Additionally in the cumulative model, the deposition of carcasses within patch stimulate increases in nitrogen, which then become available for future foraging caribou after a set leg time. If a patch is depleted of nitrogen via foraging and no nitrogen is deposited, that patch will remain depleted of nitrogen until nitrogen is deposited via a carcass or fecal matter.


## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="8OffHomCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHomCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHomCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHomCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHomCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHomCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHomCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHomCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHomCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHomCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHomCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHomCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHomCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHomCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHomCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHomCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHomCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHomCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHomCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHomCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHomCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHomCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHomCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHomCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHomCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHomCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHomCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHomCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHomCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHomCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHetCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHetCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHetCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHetCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHetCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHetCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHetCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHetCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHetCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHetCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHetCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHetCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHetCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHetCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHetCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHetCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHetCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHetCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHetCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHetCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHetCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHetCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHetCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHetCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHetCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHetCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHetCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OffHetCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OffHetCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OffHetCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHomCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHomCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHomCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHomCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHomCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHomCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHomCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHomCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHomCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHomCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHomCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHomCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHomCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHomCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHomCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHomCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHomCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHomCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHomCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHomCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHomCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHomCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHomCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHomCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHomCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHomCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHomCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHomCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHomCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHomCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHetCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHetCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHetCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHetCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHetCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHetCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHetCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHetCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHetCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHetCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHetCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHetCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHetCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHetCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHetCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHetCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHetCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHetCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHetCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHetCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHetCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHetCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHetCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHetCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHetCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHetCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHetCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="8OnHetCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "8OnHetCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "8OnHetCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHomCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHomCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHomCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHomCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHomCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHomCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHomCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHomCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHomCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHomCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHomCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHomCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHomCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHomCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHomCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHomCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHomCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHomCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHomCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHomCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHomCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHomCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHomCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHomCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHomCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHomCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHomCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHomCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHomCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHomCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHetCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHetCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHetCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHetCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHetCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHetCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHetCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHetCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHetCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHetCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHetCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHetCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHetCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHetCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHetCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHetCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHetCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHetCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHetCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHetCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHetCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHetCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHetCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHetCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHetCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHetCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHetCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OffHetCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OffHetCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OffHetCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHomCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHomCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHomCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHomCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHomCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHomCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHomCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHomCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHomCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHomCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHomCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHomCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHomCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHomCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHomCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHomCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHomCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHomCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHomCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHomCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHomCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHomCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHomCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHomCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHomCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHomCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHomCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHomCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHomCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHomCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHetCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHetCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHetCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHetCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHetCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHetCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHetCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHetCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHetCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHetCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHetCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHetCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHetCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHetCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHetCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHetCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHetCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHetCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHetCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHetCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHetCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHetCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHetCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHetCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHetCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHetCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHetCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="72OnHetCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "72OnHetCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "72OnHetCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="72"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHomCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHomCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHomCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHomCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHomCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHomCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHomCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHomCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHomCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHomCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHomCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHomCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHomCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHomCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHomCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHomCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHomCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHomCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHomCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHomCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHomCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHomCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHomCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHomCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHomCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHomCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHomCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHomCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHomCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHomCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHetCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHetCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHetCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHetCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHetCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHetCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHetCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHetCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHetCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHetCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHetCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHetCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHetCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHetCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHetCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHetCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHetCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHetCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHetCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHetCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHetCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHetCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHetCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHetCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHetCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHetCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHetCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OffHetCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OffHetCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OffHetCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHomCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHomCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHomCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHomCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHomCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHomCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHomCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHomCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHomCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHomCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHomCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHomCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHomCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHomCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHomCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHomCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHomCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHomCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHomCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHomCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHomCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHomCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHomCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHomCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHomCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHomCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHomCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHomCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHomCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHomCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHetCumulative1" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHetCumulative1Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHetCumulative1Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHetCumulative2" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHetCumulative2Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHetCumulative2Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHetCumulative3" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHetCumulative3Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHetCumulative3Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHetCumulative4" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHetCumulative4Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHetCumulative4Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHetCumulative5" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHetCumulative5Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHetCumulative5Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHetCumulative6" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHetCumulative6Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHetCumulative6Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHetCumulative7" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHetCumulative7Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHetCumulative7Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHetCumulative8" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHetCumulative8Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHetCumulative8Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHetCumulative9" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHetCumulative9Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHetCumulative9Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="180OnHetCumulative10" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup
export-world "180OnHetCumulative10Start.csv"</setup>
    <go>go</go>
    <final>export-world "180OnHetCumulative10Final.csv"</final>
    <timeLimit steps="131400"/>
    <metric>ticks</metric>
    <metric>[caribou-day] of one-of caribous</metric>
    <metric>[caribou-year] of one-of caribous</metric>
    <metric>adult-dead</metric>
    <metric>calf-dead</metric>
    <metric>N-Adult-dead-count</metric>
    <metric>N-calf-dead-count</metric>
    <metric>count caribous with [sex = "female"]</metric>
    <metric>count caribous with [sex = "male"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = true and sex = "male"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "female"]</metric>
    <metric>count caribous with [is-calf? = false and sex = "male"]</metric>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="Population">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Social-groups">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min-move-length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-exp">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Heterogeneity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Interactions">
      <value value="&quot;Cumulative&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
