################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Corr_preprocess.cpp \
../src/MainSanepicCorr_mpi.cpp \
../src/NoCorr_preprocess.cpp \
../src/map_making.cpp \
../src/parsePre.cpp \
../src/todprocess.cpp 

OBJS += \
./src/Corr_preprocess.o \
./src/MainSanepicCorr_mpi.o \
./src/NoCorr_preprocess.o \
./src/map_making.o \
./src/parsePre.o \
./src/todprocess.o 

CPP_DEPS += \
./src/Corr_preprocess.d \
./src/MainSanepicCorr_mpi.d \
./src/NoCorr_preprocess.d \
./src/map_making.d \
./src/parsePre.d \
./src/todprocess.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


