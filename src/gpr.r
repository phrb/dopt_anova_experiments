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
initial_sample <- 13
added_training_points <- 6
iterations = round((initial_budget - initial_sample) / added_training_points)
print(iterations)

repetitions <- 1000
# sd_range <- 1.96
sd_range <- 0.3
nugget <- 1e-12

results <- NULL

delta <- 0.5 # (0, 1)
scale_factor <- 30

variance_dampening <- function(dimension, delta, time) {
    return(sqrt((log((dimension * ((time * pi) ^ 2)) /
                         (6 * delta))) /
                scale_factor))
}

variance_dampening_ucb <- function(dimension, delta, time) {
    return(sqrt(log(time)))
}

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

        invisible(
            capture.output(
                gp_model <-
                    km(formula = ~ y_component_number + I(1 / y_component_number) +
                           vector_length + lws_y + I(1 / lws_y) +
                           load_overlap + temporary_size +
                           elements_number + I(1 / elements_number) +
                           threads_number + I(1 / threads_number),
                       design = select(training_sample,
                                       -row_number,
                                       -time_per_pixel),
                       response = training_sample$time_per_pixel,
                       nugget = nugget * var(training_sample$time_per_pixel),
                       #nugget = 1e-1,
                       control = list(pop.size = 400,
                                      BFGSburnin = 500))
            ))

        gp_prediction <- predict(gp_model,
                                 select(testing_sample,
                                        -row_number,
                                        -time_per_pixel),
                                 "UK")

        # testing_sample$expected_improvement <- gp_prediction$mean -
        #     (sd_range * gp_prediction$sd)

        print(c(length(training_sample[, 1]),
                variance_dampening(length(training_sample) - 2,
                                   delta,
                                   i)))

        testing_sample$expected_improvement <- gp_prediction$mean -
            (variance_dampening(length(training_sample),
                                delta,
                                length(training_sample[, 1])) * gp_prediction$sd)

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
        mutate(method = "GPR_dampening") %>%
        group_by(experiment_id)

    write.csv(best_points %>%
              filter(time_per_pixel == min(time_per_pixel)),
              "gpr_dampening_best_points.csv")

    write.csv(best_points, "gpr_dampening_all_points.csv")
}
