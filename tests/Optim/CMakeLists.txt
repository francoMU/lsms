
find_package(Armadillo REQUIRED)

#
# 1. Test
#
add_executable(optim)

target_sources(optim
        PUBLIC
        optim_test.cpp)

add_test(NAME optim_test COMMAND optim)

# Use the inline header only libary
target_include_directories(optim PUBLIC  ${OPTIM_INCLUDES})

target_include_directories(optim PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(optim PUBLIC ${ARMADILLO_LIBRARIES})

#
# 2. Test
#
add_executable(no_grad_optim)

target_sources(no_grad_optim
        PUBLIC
        no_grad_optim_test.cpp)

add_test(NAME no_grad_optim_test COMMAND no_grad_optim)

# Use the inline header only libary
target_include_directories(no_grad_optim PUBLIC  ${OPTIM_INCLUDES})

target_include_directories(no_grad_optim PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(no_grad_optim PUBLIC ${ARMADILLO_LIBRARIES})
