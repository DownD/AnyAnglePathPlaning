CFLAGS = -O3

# Required to compile if not running in HOG.
CC = g++ -DNO_HOG
EXEC = a
COMMON_CPP_FILES = main-all.cpp ScenarioLoader.cpp Timer.cpp AnyAngleAlgorithm.cpp

#MAP = 64room_000
#MAP = maze512-2-0
#MAP = random512-10-0
#MAP = AcrosstheCape
MAP = arena2
#MAP = Aurora
#MAP = Predators
#MAP = orz100d
#MAP = ost100d
#MAP = random512-40-3
#MAP = random512-40-8
#MAP = TheFrozenSea
#MAP = AR0011SR
#MAP = orz700d
#MAP = lak405d

a-euc:
	$(CC) $(CFLAGS) -DEXPERIMENT_A_EUC -o EXP_A_EUC $(COMMON_CPP_FILES) ThetaStar.cpp
	./EXP_A_EUC $(MAP)

a-oct:
	$(CC) $(CFLAGS) -DEXPERIMENT_A_OCT -o EXP_A_OCT $(COMMON_CPP_FILES) ThetaStar.cpp
	./EXP_A_OCT $(MAP)

theta:
	$(CC) $(CFLAGS) -DEXPERIMENT_T -o EXP_T $(COMMON_CPP_FILES) ThetaStar.cpp
	./EXP_T $(MAP)

lazy-theta:
	$(CC) $(CFLAGS) -DEXPERIMENT_L -o EXP_L $(COMMON_CPP_FILES) ThetaStar.cpp
	./EXP_L $(MAP)
	
theta-all:
	make a-euc
	make a-oct
	make theta
	make lazy-theta

field:
	$(CC) $(CFLAGS) -DEXPERIMENT_F -o EXP_F $(COMMON_CPP_FILES) FieldAStar.cpp
	./EXP_F $(MAP)

block:
	$(CC) $(CFLAGS) -DEXPERIMENT_B -o EXP_B $(COMMON_CPP_FILES) BlockAStar.cpp
	./EXP_B $(MAP)

sub-1-a:
	$(CC) $(CFLAGS) -DEXPERIMENT_SUB_1_A -o EXP_SUB_1_A $(COMMON_CPP_FILES) SubgoalAA.cpp
	./EXP_SUB_1_A $(MAP)

sub-1-t:
	$(CC) $(CFLAGS) -DEXPERIMENT_SUB_1_T -o EXP_SUB_1_T $(COMMON_CPP_FILES) SubgoalAA.cpp
	./EXP_SUB_1_T $(MAP)	

sub-2-a:
	$(CC) $(CFLAGS) -DEXPERIMENT_SUB_2_A -o EXP_SUB_2_A $(COMMON_CPP_FILES) SubgoalAA.cpp
	./EXP_SUB_2_A $(MAP)

sub-2-t:
	$(CC) $(CFLAGS) -DEXPERIMENT_SUB_2_T -o EXP_SUB_2_T $(COMMON_CPP_FILES) SubgoalAA.cpp
	./EXP_SUB_2_T $(MAP)

sub-n-a:
	$(CC) $(CFLAGS) -DEXPERIMENT_SUB_10000_A -o EXP_SUB_10000_A $(COMMON_CPP_FILES) SubgoalAA.cpp
	./EXP_SUB_10000_A $(MAP)

sub-n-t:
	$(CC) $(CFLAGS) -DEXPERIMENT_SUB_10000_T -o EXP_SUB_10000_T $(COMMON_CPP_FILES) SubgoalAA.cpp
	./EXP_SUB_10000_T $(MAP)

sub-all:
	make sub-1-a
	make sub-1-t
	make sub-2-a
	make sub-2-t
	make sub-n-a
	make sub-n-t

anya:
	$(CC) $(CFLAGS) -DEXPERIMENT_ANYA -o EXP_ANYA $(COMMON_CPP_FILES) ANYA.cpp
	./EXP_ANYA $(MAP)

all:
	make a-euc
	make a-oct
	make theta
	make lazy-theta
	make field
	make block
	make sub-1-a
	make sub-1-t
	make sub-2-a
	make sub-2-t
	make sub-n-a
	make sub-n-t
	make anya
	make move-data

test-all:
	./EXP_A_EUC $(MAP)
	./EXP_A_OCT $(MAP)
	./EXP_T $(MAP)
	./EXP_L $(MAP)
	./EXP_F $(MAP)
	./EXP_B $(MAP)
	./EXP_SUB_1_A $(MAP)
	./EXP_SUB_1_T $(MAP)	
	./EXP_SUB_2_A $(MAP)
	./EXP_SUB_2_T $(MAP)
	./EXP_SUB_10000_A $(MAP)
	./EXP_SUB_10000_T $(MAP)
	./EXP_ANYA $(MAP)
	make move-data

move-data:
	mkdir -p $(MAP)
	mv *-$(MAP) $(MAP)

