# Specify which compiler to use:
CXX=g++

# Use CFLAGS to specify additional compiler options to use:
C11FLAGS= -std=c++11

# Define variables:
EXE_FILE= exekmc
SOURCE_FILES= kmc_*.cpp 
ALL_FILES= kmc_*

# Build executable:
$(EXE_FILE) : $(ALL_FILES)
	$(CXX) -o $(EXE_FILE) $(C11FLAGS) $(SOURCE_FILES)
	chmod +x $(EXE_FILE)	# turn on execute permissions

# Rules:
clean:
	rm $(EXE_FILE)

