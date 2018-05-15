library(AlgDesign)
library(dplyr)

check_prediction <- function(regression, prediction, target_size) {
  predicted_effects <- summary(regression)[[1]][["Pr(>F)"]]
  names(predicted_effects) <- trimws(rownames(summary(regression)[[1]]))
  predicted_effects <- data.frame(t(predicted_effects), check.names = FALSE)

  ordered_effects <- names(predicted_effects[ , order(predicted_effects[1, ])])
  return(prediction %in% ordered_effects[1:target_size])
}

complete_data = read.csv("../data/search_space.csv", header = TRUE)
# str(complete_data)

accuracies_file <- "../data/testing_lm_bm_accuracies.csv"
read_file <- "../data/old_complete_1000.csv"
write_file <- "../data/testing_1000.csv"

results <- read.csv(read_file, strip.white=T, header=T)

budget <- 120

iterations <- 1000

factors = c("elements_number", "y_component_number",
            "vector_length", "temporary_size",
            "load_overlap", "threads_number",
            "lws_y")

accuracies <- data.frame(first_step = 0, second_step = 0, third_step = 0)

for (i in 1:iterations) {
    str(i)

    used <- 0

    sample_size <- 50
    data <- complete_data[, c(factors, "time_per_pixel")]
    design <- data[sample(rownames(data), sample_size, replace = FALSE), ]
    #str(design)

    experiments <- rownames(design)

    # Complete model:
    regression <- aov(time_per_pixel ~ y_component_number + I(1 / y_component_number) +
                                       vector_length + lws_y + I(1 / lws_y) +
                                       load_overlap + temporary_size +
                                       elements_number + I(1 / elements_number) +
                                       threads_number + I(1 / threads_number),
                      data = design)

    if (all(check_prediction(regression, c("vector_length", "lws_y"), 3)) |
        all(check_prediction(regression, c("vector_length", "I(1/lws_y)"), 3))) {
        accuracies$first_step <- accuracies$first_step + 1
    }

    used <- used + sample_size

    #print(summary(regression))

    # Checking the ANOVA summary we can identify at least two variables
    # that seem to have greater impact: 'vector_length' and 'lws_y'.
    # Let's fix those variables to their best predicted value so far,
    # then fit a new model without them

    predicted_best <- data[predict(regression, data) == min(predict(regression, data)), ]
    # predicted_best

    data <- complete_data[complete_data$vector_length == predicted_best$vector_length &
                          complete_data$lws_y == predicted_best$lws_y, c(factors, "time_per_pixel")]

    sample_size <- 40
    unique_experiments <- data[rownames(data)[!(rownames(data) %in% experiments)], ]
    design <- unique_experiments[sample(rownames(unique_experiments), sample_size, replace = FALSE), ]
    #str(design)

    experiments <- c(experiments, rownames(design))

    used <- used + sample_size

    regression <- aov(time_per_pixel ~ y_component_number + I(1 / y_component_number) +
                                       load_overlap + temporary_size +
                                       elements_number + I(1 / elements_number) +
                                       threads_number + I(1 / threads_number),
                      data = design)

    if (all(check_prediction(regression, c("y_component_number", "threads_number"), 4)) |
        all(check_prediction(regression, c("I(1/y_component_number)", "threads_number"), 4)) |
        all(check_prediction(regression, c("y_component_number", "I(1/threads_number)"), 4)) |
        all(check_prediction(regression, c("I(1/y_component_number)", "I(1/threads_number)"), 4))) {
        accuracies$second_step <- accuracies$second_step + 1
    }

    #print(summary(regression))

    # Checking the ANOVA summary we can identify at least two variables
    # that seem to have greater impact: 'y_component_number' and 'threads_number'.
    # Let's fix those variables to their best predicted value so far,
    # then fit a new model without them

    predicted_best <- data[predict(regression, data) == min(predict(regression, data)), ]
    # predicted_best

    data <- complete_data[complete_data$vector_length == predicted_best$vector_length &
                          complete_data$lws_y == predicted_best$lws_y &
                          complete_data$y_component_number == predicted_best$y_component_number &
                          complete_data$threads_number == predicted_best$threads_number, c(factors, "time_per_pixel")]

    sample_size <- 20
    unique_experiments <- data[rownames(data)[!(rownames(data) %in% experiments)], ]

    if (nrow(unique_experiments) > sample_size) {
        design <- unique_experiments[sample(rownames(unique_experiments), sample_size, replace = FALSE), ]
    } else {
        design <- unique_experiments
    }

    #str(design)
    experiments <- c(experiments, rownames(design))

    used <- used + nrow(design)

    regression <- aov(time_per_pixel ~ load_overlap + temporary_size +
                        elements_number + I(1 / elements_number),
                      data = design)

    if (all(check_prediction(regression, c("elements_number"), 2)) |
        all(check_prediction(regression, c("I(1/elements_number)"), 2))) {
        accuracies$third_step <- accuracies$third_step + 1
    }

    #print(summary(regression))

    # Checking the ANOVA summary we can identify, at last, one variable
    # that seem to have greater impact: 'elements_number'
    # and 'elements_number'.
    # Let's fix it to their best predicted value so far,
    # then fit a new model without it

    predicted_best <- data[predict(regression, data) == min(predict(regression, data)), ]
    # predicted_best

    data <- complete_data[complete_data$vector_length == predicted_best$vector_length &
                          complete_data$lws_y == predicted_best$lws_y &
                          complete_data$y_component_number == predicted_best$y_component_number &
                          complete_data$threads_number == predicted_best$threads_number &
                          complete_data$elements_number == predicted_best$elements_number, c(factors, "time_per_pixel")]

    sample_size <- 5
    unique_experiments <- data[rownames(data)[!(rownames(data) %in% experiments)], ]

    if (nrow(unique_experiments) > sample_size) {
        design <- unique_experiments[sample(rownames(unique_experiments), sample_size, replace = FALSE), ]
    } else {
        design <- unique_experiments
    }

    if (nrow(design) > 0) {
        str(design)
        experiments <- c(experiments, rownames(design))

        used <- used + nrow(design)

        regression <- aov(time_per_pixel ~ load_overlap + temporary_size,
                          data = design)

        # print(summary(regression))

        predicted_best <- data[predict(regression, data) == min(predict(regression, data)), ]
    }

    best <- complete_data[complete_data$time_per_pixel == min(complete_data$time_per_pixel), ]
    best_row <- rownames(best)

    predicted_best$slowdown <- predicted_best$time_per_pixel / best$time_per_pixel
    predicted_best$method <- rep("LM_bm", nrow(predicted_best))
    predicted_best$point_number <- rep(used, nrow(predicted_best))
    predicted_best$vector_recompute <- rep("true", nrow(predicted_best))

    predicted_best <- predicted_best[, c("elements_number", "y_component_number",
                                        "vector_length", "temporary_size", "vector_recompute",
                                        "load_overlap", "threads_number", "lws_y",
                                        "time_per_pixel", "point_number", "method",
                                        "slowdown")]
    #print(predicted_best)
    results <- rbind(results, predicted_best)
}

accuracies <- accuracies / iterations

write.csv(accuracies, accuracies_file, row.names = FALSE)
write.csv(results, write_file, row.names = FALSE)
