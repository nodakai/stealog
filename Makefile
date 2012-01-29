all : grad

#CXXFLAGS= -Wall -Wextra -O0 -g3
CXXFLAGS= -Wall -Wextra -Ofast -g3

grad : grad.cpp Makefile
	$(CXX) $< -o $@ $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS)

merged.txt : result2.txt merge.rb
	ruby merge.rb $< > $@
