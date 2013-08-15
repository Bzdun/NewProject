﻿CREATE TABLE [dbo].[Staging_ExhibitView] (
    [ID]                UNIQUEIDENTIFIER   NOT NULL,
    [Title]             NVARCHAR (100)     NOT NULL,
    [ThresholdID]       UNIQUEIDENTIFIER   NOT NULL,
    [RegimeID]          UNIQUEIDENTIFIER   NOT NULL,
    [ContentDate]       DATETIMEOFFSET (7) NULL,
    [ContentYear]       DECIMAL (8, 4)     NULL,
    [CreatedOn]         SMALLDATETIME      NULL,
    [CreatedBy]         UNIQUEIDENTIFIER   NULL,
    [ModifiedOn]        SMALLDATETIME      NULL,
    [ModifiedBy]        UNIQUEIDENTIFIER   NULL,
    [IsVisible]         BIT                NULL,
    [IsDeleted]         BIT                NULL,
    [TimeUnitID]        UNIQUEIDENTIFIER   NULL,
    [UniqueID]          INT                NOT NULL,
    [Sequence]          INT                NULL,
    [CurrVersion]       INT                NULL,
    [Version_ExhibitID] UNIQUEIDENTIFIER   NOT NULL,
    [IncludeInVersion]  INT                NULL,
    [ExhibitViewID]     UNIQUEIDENTIFIER   NOT NULL,
    PRIMARY KEY CLUSTERED ([ExhibitViewID] ASC)
);

