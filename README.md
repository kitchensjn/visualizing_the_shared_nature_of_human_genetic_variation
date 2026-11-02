# Visualizing the shared nature of human genetic variation
This folder contains all of the code to recreate the figures for the manuscript Visualizing the shared nature of human genetic variation by Kitchens, Angel-Clark, and Coop (2026).

## Online Activity | [Demo](https://james-kitchens.com/visualizing_the_shared_nature_of_human_genetic_variation/)

We've developed a short online assignment to introduce students to Euler diagrams and their application to genetic data. If you would like to record your students' responses as part of a larger class discussion, follow the below instructions:

1) If you haven't yet, fork this repository on GitHub
2) Set up a Google Sheet and Google Apps Script to receive your students' responses:
    - [General instructions](https://github.com/jamiewilson/form-to-google-sheets) for set up
    - Add these column names to Row 1: `timestamp`, `time_spent`, `finished_codeword`, `interpretability_overlap0`, `interpretability_overlap1`, `build_your_own_overlap0`, `build_your_own_overlap1`, `thousand_genomes_overlap0`, `thousand_genomes_overlap1`, `thousand_genomes_overlap2`
    - Paste the link to your Google Sheet in Line 575 of `docs/_layouts/home.html`
    - Commit these changes to your GitHub fork
3) This activity is created with Jekyll and can be hosted through GitHub Pages:
    - Within the GitHub repository, go to `Settings` and then `Pages` in the left menu bar
    - In Build and deployment, make sure that `Source` is set to `Deploy from a branch` and `Branch` is set to `docs/`
    - In a minute, the activity should deploy. The link will be at the top of your current page.
    - Make sure to test that the answers are being recorded to your Google Sheet!

If you have any questions, feel free to reach contact me (**James Kitchens**) at kitchensjn at gmail.com.

