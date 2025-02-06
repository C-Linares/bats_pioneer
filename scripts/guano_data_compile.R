# Load required libraries
library(future)
library(future.apply)
library(guano)

# Set up a plan for asynchronous processing
plan(multisession)

# Define the function to be run in the background
process_files <- function() {
  execution_time <- system.time({
    dt2 <- read.guano.dir(dirname = 'Z:/PioneerLights_2021/', recursive = T)
  })
  return(execution_time)
}

# Run the function as a background job
future_result <- future({
  process_files()
})

# Print the execution time once the job is complete
execution_time <- value(future_result)
print(execution_time)


