Feature: Change is always occurring in the Earth system
  The workbench must say what is changing, where, by how much, and whether the
  change is statistically significant.

  Scenario: A warming series is reported as significant change
    Given a monthly series warming at 0.03 per month
    When the workbench analyses it
    Then a positive trend slope is reported
    And the trend is declared significant
    And the trend estimator ensemble reports a primary method

  Scenario: A flat series is not declared significant
    Given a monthly series with no trend
    When the workbench analyses it
    Then the trend is not declared significant

  Scenario: The spatial field localises where change happens
    When the spatial grid is requested
    Then every cell has a trend slope
    And at least one region is flagged as a hotspot

  Scenario: The formula lab runs every selected formula at once
    Given the formula lab catalogue
    When a batch of three formulas is run
    Then three results are returned
    And each result either has a value or a documented error

  Scenario: Live NASA data can be pulled with provenance
    When the GISTEMP series is pulled
    Then more than a thousand samples are returned
    And the provenance records the NASA endpoint
