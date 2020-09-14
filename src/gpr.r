library(DiceKriging)
library(dplyr)

complete_data <- read.csv("dopt_anova_experiments/data/search_space.csv",
                          header = TRUE)
complete_data <- complete_data %>%
    select(-vector_recompute) %>%
    mutate(row_number = row_number(),
           load_overlap = as.numeric(as.factor(load_overlap)))

str(complete_data)

global_optimum <- filter(complete_data, time_per_pixel == min(time_per_pixel))

initial_budget <- 120
initial_sample <- 24
added_training_points <- 6
iterations = (initial_budget - initial_sample) / added_training_points

repetitions <- 1000

results <- NULL

for(j in 1:repetitions){
    testing_sample <- complete_data
    training_sample <- NULL

    for(i in 1:iterations){
        if(is.null(training_sample)){
            training_sample <- slice_sample(testing_sample,
                                            n = initial_sample)
        }

        testing_sample <- testing_sample %>%
            filter(!(row_number %in% training_sample$row_number))

        gp_model <- km(formula = ~ .,
                       design = select(training_sample,
                                       -row_number,
                                       -time_per_pixel),
                       response = training_sample$time_per_pixel,
                       nugget = 1e-8 * var(training_sample$time_per_pixel),
                       control = list(pop.size = 400,
                                      BFGSburnin = 500))

        gp_prediction <- predict(gp_model,
                                 select(testing_sample,
                                        -row_number,
                                        -time_per_pixel),
                                 "UK")

        testing_sample$expected_improvement <- gp_prediction$mean - (1.96 * gp_prediction$sd)

        new_training_sample <- testing_sample %>%
            arrange(expected_improvement)

        testing_sample <- select(testing_sample, -expected_improvement)

        new_training_sample <- select(new_training_sample[1:added_training_points, ],
                                      -expected_improvement)

        training_sample <- bind_rows(training_sample,
                                     new_training_sample)
    }

    training_sample <- training_sample %>%
        mutate(measurement_order = row_number(),
               experiment_id = j,
               slowdown = time_per_pixel /
                   global_optimum$time_per_pixel)

    if(is.null(results)){
        results <- training_sample
    } else{
        results <- bind_rows(results,
                             training_sample)
    }

    best_points <- results %>%
        mutate(method = "GPR") %>%
        group_by(experiment_id)

    write.csv(best_points %>%
              filter(time_per_pixel == min(time_per_pixel)),
              "best_points.csv")

    write.csv(best_points, "all_points.csv")
}
