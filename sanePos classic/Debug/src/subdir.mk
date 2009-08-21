################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/MainSanepicCorr_mpi.cpp \
../src/blastSpecific.cpp \
../src/parsePos.cpp \
../src/sanePos_Preprocess.cpp \
../src/sanePos_mapmaking.cpp 

OBJS += \
./src/MainSanepicCorr_mpi.o \
./src/blastSpecific.o \
./src/parsePos.o \
./src/sanePos_Preprocess.o \
./src/sanePos_mapmaking.o 

CPP_DEPS += \
./src/MainSanepicCorr_mpi.d \
./src/blastSpecific.d \
./src/parsePos.d \
./src/sanePos_Preprocess.d \
./src/sanePos_mapmaking.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -UMPICH_IGNORE_CXX_SEEK -I"/home/matthieu/workspace/nr/src" -I"/home/matthieu/workspace/Sanelib/src" -I"/home/matthieu/workspace/saneIO/src" -I"/home/matthieu/workspace/iniparser/src" -O3 -p -pg -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


