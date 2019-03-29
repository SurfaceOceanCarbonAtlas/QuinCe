CREATE TABLE sensor_values (
  id BIGINT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  dataset_id INT NOT NULL,
  file_column INT NOT NULL,
  date BIGINT(20) NOT NULL,
  value VARCHAR(100) NULL,
  auto_qc_flag SMALLINT(2) DEFAULT -1000,
  auto_qc_message VARCHAR(255),
  user_qc_flag SMALLINT(2) DEFAULT -1000,
  user_qc_message VARCHAR(255),
  CONSTRAINT SENSORVALUE_DATASET
    FOREIGN KEY (dataset_id)
    REFERENCES dataset (id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION);

    