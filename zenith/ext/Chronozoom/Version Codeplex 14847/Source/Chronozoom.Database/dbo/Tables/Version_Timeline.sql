﻿CREATE TABLE [dbo].[Version_Timeline] (
    [ID]                 UNIQUEIDENTIFIER   NOT NULL,
    [Title]              NVARCHAR (100)     NOT NULL,
    [ThresholdID]        UNIQUEIDENTIFIER   NOT NULL,
    [RegimeID]           UNIQUEIDENTIFIER   NOT NULL,
    [FromContentDate]    DATETIMEOFFSET (7) NULL,
    [FromContentYear]    DECIMAL (8, 4)     NULL,
    [ToContentDate]      DATETIMEOFFSET (7) NULL,
    [ToContentYear]      DECIMAL (8, 4)     NULL,
    [ParentTimelineID]   UNIQUEIDENTIFIER   NULL,
    [CreatedOn]          SMALLDATETIME      NULL,
    [CreatedBy]          UNIQUEIDENTIFIER   NULL,
    [ModifiedOn]         SMALLDATETIME      NULL,
    [ModifiedBy]         UNIQUEIDENTIFIER   NULL,
    [IsVisible]          BIT                NULL,
    [IsDeleted]          BIT                NULL,
    [FromTimeUnit]       UNIQUEIDENTIFIER   NULL,
    [ToTimeUnit]         UNIQUEIDENTIFIER   NULL,
    [SourceURL]          NVARCHAR (500)     NULL,
    [Attribution]        NVARCHAR (500)     NULL,
    [UniqueID]           INT                NOT NULL,
    [Sequence]           INT                NULL,
    [Height]             DECIMAL (8, 4)     NULL,
    [CurrVersion]        INT                NULL,
    [Version_TimelineID] UNIQUEIDENTIFIER   NOT NULL,
    PRIMARY KEY CLUSTERED ([Version_TimelineID] ASC)
);

